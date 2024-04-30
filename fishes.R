library(here)
library(fs)
library(tidyverse)
library(janitor)
library(RColorBrewer)
library(patchwork)
library(vegan)
library(paletteer)
library(httr)
library(EcolUtils)
library(dendextend)
library(ggtext)


# setup -------------------------------------------------------------------

# whether to save to pdf all the time
save_pdf <- FALSE

# source utility functions
source("util.R")

# make figures directory
fig_dir <- here("output","figures")
dir_create(fig_dir,recurse = TRUE)

# get rid of annoying '`summarise()` has grouped output by' message
options(dplyr.summarise.inform = FALSE)

# get base data directory
data_dir <- here("data")

# read markers
markers <- read_lines(path(data_dir,"markers.txt"))

# make a consistent color palette for markers
marker_pal <- paletteer_d("ggthemes::Miller_Stone",n=length(markers)) %>%
  as.character() %>%
  set_names(markers)

# fwd/reverse maps of cleaned marker names to their display-worthy versions
# (useful for plot labels and so forth)
mmap <- markers %>%
  set_names(make_clean_names(markers))
mmap_r <- markers %>%
  make_clean_names() %>%
  set_names(markers)

# load metadata
metadata <- read_csv(path(data_dir,"metadata.csv")) %>%
  clean_names() %>%
  mutate(clean_sample = make_clean_names(sample))

# map of regex patterns to match the various sample types
dataset_map <- list(
  sharkpen = "_sp[0-9]+$",
  aquarium = "_wa[0-9]+$",
  mock = "_mc[0-9]+$"
)

# map dataset names to display titles
title_map <- c(
  sharkpen = "Shark pen",
  aquarium = "Waikīkī Aquarium",
  mock = "Mock community"
)

# specific taxonomic orders to filter out
filter_orders <- c("Artiodactyla", "Carnivora", "Mammalia", "Phlebobranchia", "Stolidobranchia", "Primates")

# minimum total sample abundance
min_total <- 25

# minimum read count
min_reads <- 10

# number of rarefaction permutations
rarefy_perm <- 100

# whether to filter out unidentified taxa
filter_unid <- TRUE

# function to produce summary statistics for plotting
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  
  ymin <- ifelse(is.na(ymin),m,ymin)
  ymax <- ifelse(is.na(ymax),m,ymax)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# shortcut to do sum(x,na.rm=TRUE)
nasum <- function(x) sum(x,na.rm=TRUE)


# download ncbi taxonomy
get_lineage <- function(nlf) {
  # nlf <- path(data_dir,"rankedlineage.dmp")
  if (!file_exists(nlf)) {
    url <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"
    taxdump <- path(data_dir,"taxdump.zip")
    resp <- GET(url,write_disk(taxdump,overwrite = T),progress())
    if (resp$status_code == 200) {
      unzip(taxdump,files="rankedlineage.dmp",exdir=data_dir)       
    } else {
      return(NULL)
    }
  }
  read_tsv(nlf,col_types = "i_c_c_c_c_c_c_c_c_c_",
           col_names = c("taxid","taxon","species","genus","family","order","class","phylum","kingdom","domain")) %>%
    select(taxid,taxon,domain,kingdom,phylum,class,order,family,genus,species,taxon)
}

# coalesce across
coacross <- function(...) {
  coalesce(!!!across(...))
}

# initial data import & filtering -----------------------------------------

# all sample column names
all_samples <- c()

# load raw data and do some pre-filtering
# map through markers and return a concatenated tibble
fishes <- markers %>%
  map_dfr(~{
    # save marker name into something nicer than `.x`
    marker <- .x
    
    # read marker data and filter to just chordates
    chordates <- read_tsv(path(data_dir,str_glue("{marker}_data.tsv")),col_types = cols()) %>%
      # make nicer 'machine-readable' column names
      clean_names() %>%
      filter(phylum == "Chordata")
    
    # subtract reads in blanks and other stuff
    noblanks <- chordates %>%
      mutate(
        # replace "dropped" with blank string
        across(class:species,~replace(.x,which(.x == "dropped"),"")),
        # sum across blanks into column `blanks`
        blanks = rowSums(pick(matches("blank")),na.rm=TRUE),
        # subtract blank reads from sample reads
        across(-c(domain:seq_length,blanks),~.x-blanks),
        # set negative reads to zero
        across(where(is.numeric),~replace(.x,which(.x < 0),0))
      ) %>%
      # get rid of blank and sum columns (anything with the text "blank" in it)
      select(-matches("blank"))
    
    # pull out sample IDs
    samples <- noblanks %>%
      select(-c(domain:seq_length)) %>%
      names()
    
    # concatenate sample names to all samples (<<- does so in the global context)
    all_samples <<- c(all_samples,samples)
    
    # make ID column name
    id_col <- make_clean_names(str_c(marker,"_id"))
    
    # finally, collapse by species and return result
    noblanks %>% 
      group_by(across(domain:species)) %>%
      # group_by(domain,kingdom,phylum,class,order,family,genus,species) %>%
      summarise(
        # sum sample reads
        across(all_of(samples),nasum),
        # pull out representative zotu (first one), and specify which marker it is
        representative = str_glue("{marker}:{otu[1]}"),
        
        # concatenate zotus with marker name like this marker(zotu1,zotu2,zotu3,...)
        zotus = str_glue("{marker}({str_c(otu,collapse=',')})")
      ) %>%
      ungroup() %>%
      select(domain:species,representative,zotus,all_of(samples))
  }) %>%
  # these operations are done on the fully concatenated dataset
  mutate(
    # replace NA taxa with blank string
    across(domain:species,~replace_na(.x,"")),
    # replace NA reads with zeroes
    across(where(is.numeric),~replace_na(.x,0))
  ) %>%
  # sum reads by taxon and concatenate zotus/representatives
  group_by(across(domain:species)) %>%
  summarise(
    # concatenate representative zotus
    representative=str_c(representative,collapse=','),
    # concatenate zotu lists
    zotus=str_c(zotus,collapse=','),
    # sum read counts ()
    across(all_of(all_samples),nasum)
  ) %>%
  ungroup()

# pull out sample names
samples <- all_samples
# samples <- fishes %>%
#   select(-c(domain:species)) %>%
#   names()


# some final data wrangling and filtering
fishes_filtered <- fishes %>%
  # create an ID column (we call it 'zotu', but it's not really a zotu)
  mutate(zotu = str_glue("taxon{row_number()}")) %>%
  # properly order columns
  select(domain:species,zotu,everything()) %>%
  # filter out undesired taxonomic orders
  filter(!order %in% filter_orders) %>%
  # replace reads below minimum threshold with zero
  mutate(across(all_of(samples),~replace(.x,which(.x < min_reads),0))) %>%
  # get rid of samples with total read counts over minimum threshold
  select(c(domain:species,zotu,representative,zotus),where(~is.numeric(.x) && sum(.x) >= min_total))



# use the map to break out the different sample types and pivot to long format
# do one that's rarefied and one that's not
datasets <- c(FALSE,TRUE) %>%
  set_names(c("raw","rarefied")) %>%
  map(~{
    rarefy <- .x
    dataset_map %>%
      map(~{
        ff <- fishes_filtered %>%
          # grab samples that match the current sample type
          select(domain:species,zotu,matches(.x))
        # rarefy datasets if so desired
        if (rarefy) {
          # pull out taxonomy part to stick back on later
          taxonomy <- ff %>% 
            select(zotu,domain:species) 
          ff <- ff %>%
            # get only zotu and samples
            select(-c(domain:species)) %>%
            # swap to zotu x samples wide format
            pivot_longer(-zotu,names_to="sample",values_to="reads") %>%
            pivot_wider(names_from="zotu",values_from="reads",values_fill=0) %>%
            # set sample as rownames
            column_to_rownames("sample") %>%
            # do rarefaction
            rrarefy.perm(n=rarefy_perm) %>%
            # convert back to tibble, make row names the `sample` column
            as_tibble(rownames="sample") %>%
            # convert back to sample x zotu wide format
            pivot_longer(-sample,names_to="zotu",values_to="reads") %>%
            pivot_wider(names_from="sample",values_from="reads") %>%
            # join taxonomy back
            left_join(taxonomy,by="zotu") %>%
            # order columns correctly
            select(domain:species,zotu,everything())
        }
        ff %>%
          # switch to long format
          pivot_longer(-c(domain:species,zotu),names_to="sample",values_to="reads") %>%
          # join in sample metadata
          left_join(metadata,by=c("sample" = "clean_sample"),suffix=c("","_display")) %>%
          # group by zotu (which I guess it probably is by default)
          group_by(domain,kingdom,phylum,class,order,family,genus,species,marker,type,sample,zotu) %>%
          # and sum up reads (we probably don't actually need to do this)
          summarise(reads = sum(reads)) %>%
          ungroup() %>%
          # calculate relative read abundance for each marker
          group_by(marker) %>%
          mutate(rel = reads/sum(reads)) %>%
          ungroup() %>% 
          # filter out zero reads
          filter(reads > 0) %>%
          # set blank taxa to unidentified
          mutate(
            across(domain:species,~replace(.x,which(.x == ""),"unidentified")),
            marker = factor(marker,levels=markers)
          )
      })
  })

# start here

# taxon bar plots ---------------------------------------------------------

# taxonomic levels at which to group bar plots
plot_levels <- c("family","order","class")

# create a list of bar plots by plot level and dataset
rel_taxon_plotz <- datasets %>%
  map(~{
    dataset <- .x
    #map through plot levels
    plot_levels %>%
      set_names() %>%
      map(~{
        # save plot level into `pl`
        pl <- .x
        # make consistent taxon palette across datasets
        # pull all possible values for the plot level
        groups <- dataset %>%
          map(~.x %>% pull(all_of(pl))) %>%
          unlist() %>%
          sort() %>%
          unique()    
        # make palette generator and generate a palette
        pal = colorRampPalette(brewer.pal(12, "Paired"))(length(groups))
        # name the palette
        names(pal) <- groups
        # map through datasets
        dataset %>%
          imap(~{
            # save dataset name
            nn <- .y
            # sum up relative reads for marker/plot level combos
            dd <- .x %>%
              filter(!filter_unid | .data[[pl]] != "unidentified") %>%
              # recalculate rel in case we've filtered out unidentified things
              group_by(marker) %>%
              mutate(rel = reads/sum(reads)) %>%
              ungroup() %>% 
              # sum rel for marker and plot level
              group_by(marker,across(all_of(pl))) %>%
              summarise(rel = sum(rel)) %>%
              ungroup() 
            # do plotting
            ggplot(dd) + 
              geom_col(aes(x=marker,y=rel,fill=.data[[pl]])) +
              scale_x_discrete(drop = TRUE) +
              scale_y_continuous(expand = c(0, 0), limits = c(0,1.01) )+
              scale_fill_manual(values=pal,drop=FALSE,name=str_to_sentence(pl)) + 
              theme_bw() +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    axis.text=element_text(size=12),
                    panel.border = element_blank(),
                    axis.title.y = element_text(size=14),
                    axis.line.x = element_line(color="black", linewidth = .5),
                    axis.line.y = element_line(color="black", linewidth = .5)) +
              labs(title=title_map[nn],y=str_glue("Relative {pl} abundance"),x="Marker")
          })
      }) 
  })



# reduce all these plots to one big one using patchwork
# (it looks like shit because there's too much going on)

# make composite plots for each plot level
rel_taxon_composites <- rel_taxon_plotz %>%
  map(~{
    .x %>%
      imap(~{
        .x %>%
          map(~.x + theme(axis.text.x = element_text(angle=45,hjust=1))) %>%
          reduce(`+`) +
          plot_annotation(title=.y)
      })
  })

# see them like this
# rel_taxon_composites$raw$family
# rel_taxon_composites$raw$order

# save them
if (save_pdf) {
  rel_taxon_plotz %>%
    iwalk(~{
      r <- .y
      .x %>% iwalk(~{
        pl <- .y
        .x %>% iwalk(~{
          p <- .y
          ggsave(path(fig_dir,str_glue("taxon_{p}_{pl}_{r}.pdf")),.x,device=cairo_pdf,width=14,height=12,units="in")
        })
      })
    })
}

# zotu bar plots ----------------------------------------------------------

# plot zotus by taxonomic level
plot_levels <- c("family","species")

# are zotu plots relative?
zotu_rel <- FALSE

# create a list of bar plots by plot level and dataset
# map through plot levels
zotu_plotz <- datasets %>%
  map(~{
    dataset <- .x
    plot_levels %>%
      set_names() %>%
      map(~{
        # save plot level
        pl <- .x
        # map through datasets
        # make consistent taxon palette across datasets
        # pull all possible values for the plot level
        groups <- dataset %>%
          map(~.x %>% pull(all_of(pl))) %>%
          unlist() %>%
          sort() %>%
          unique()    
        # make palette generator and generate a palette
        pal = colorRampPalette(brewer.pal(12, "Paired"))(length(groups))
        # name the palette
        names(pal) <- groups
        dataset %>%
          imap(~{
            nn <- .y
            # calculate relative zotu abundance (not read abundance)
            # for marker and plot level
            dd <- .x %>%
              filter(!filter_unid | .data[[pl]] != "unidentified") %>%
              count(marker,across(all_of(pl))) %>%
              group_by(marker) %>%
              mutate(rel = n / sum(n)) %>%
              ungroup()
            if (!zotu_rel)
              dd <- dd %>% mutate(rel=n)
            # do plotting
            ggplot(dd) + 
              geom_col(aes(x=marker,y=rel,fill=.data[[pl]])) + 
              scale_fill_manual(values=pal,name=str_to_sentence(pl)) + 
              theme_bw() +
              theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text = element_text(size = 12),
                panel.border = element_blank(),
                axis.title.y = element_text(size = 14),
                axis.line.x = element_line(color = "black", linewidth = 0.5),
                axis.line.y = element_line(color = "black", linewidth = 0.5),
                legend.key.size = unit(1.2, "lines"),  # Adjust the legend key size here
                legend.text = element_text(size = 14)   # Adjust the legend text size here
              ) +
              labs(title=title_map[nn],y="zOTUs",x="Marker")
          })
      })
  })

# make composite plots for each plot level, but they'll be hard to read
zotu_composites <- zotu_plotz %>%
  map(~{
    .x %>%
      imap(~{
        .x %>%
          map(~.x + theme(axis.text.x = element_text(angle=45,hjust=1))) %>%
          reduce(`+`) +
          plot_annotation(title=.y)
      })
  })

# the family-level composite is ok, but the species one is overwhelmed by legend so we'll hide it
# zotu_composites$raw$family
# zotu_composites$raw$species & theme(legend.position = "none")

# save them
if (save_pdf) {
  zotu_plotz %>%
    iwalk(~{
      r <- .y
      .x %>% iwalk(~{
        pl <- .y
        .x %>% iwalk(~{
          p <- .y
          ggsave(path(fig_dir,str_glue("zotu_{p}_{pl}_{r}.pdf")),.x,device=cairo_pdf,width=14,height=12,units="in")
        })
      })
    })
}


# relative zotu abundance heatmaps ----------------------------------------

# plot zotus by taxonomic level
plot_levels <- c("family","species")

rel_zotu_plotz <- datasets %>%
  map(~{
    dataset <- .x
    plot_levels %>%
      set_names() %>%
      map(~{
        pl <- .x
        dataset %>%
          imap(~{
            dsn <- .y
            dd <- .x %>%
              filter(!filter_unid | .data[[pl]] != "unidentified") %>%
              # recalculate rel in case we've filtered out unidentified things
              group_by(marker) %>%
              mutate(rel = reads/sum(reads)) %>%
              ungroup() %>% 
              # sum rel for marker and plot level
              group_by(marker,across(all_of(pl))) %>%
              summarise(rel = sum(rel), reads = sum(reads)) %>%
              ungroup() %>%
              pivot_wider(id_cols=-reads,names_from = "marker", values_from = "rel", values_fill = 0) %>%
              pivot_longer(any_of(markers),names_to = "marker", values_to = "rel") %>%
              # mutate("{pl}" := fct_relevel(.data[[pl]],"unidentified",after=Inf)) %>%
              mutate("{pl}" := fct_relevel(.data[[pl]],"unidentified")) %>%
              suppressWarnings()
            
            ggplot(dd) + 
              geom_tile(aes(x=marker,y=.data[[pl]],fill=rel),color="black") + 
              scale_fill_paletteer_c("viridis::cividis")  +
              scale_y_discrete(limits=rev) + 
              labs(title=title_map[[dsn]],x="Marker",y=str_to_sentence(pl))
          })
      })
  })

# save to pdfs
if (save_pdf) {
  rel_zotu_plotz %>%
    iwalk(~{
      r <- .y
      .x %>% iwalk(~{
        pl <- .y
        .x %>% iwalk(~{
          p <- .y
          ggsave(path(fig_dir,str_glue("rel_zotu_{p}_{pl}_{r}.pdf")),.x,device=cairo_pdf,width=14,height=12,units="in")
        })
      })
    })
}

# box plots for diversity metrics -----------------------------------------

# map through datasets
div_plotz <- datasets %>%
  map(~{
    dataset <- .x
    .x %>%
      imap(~{
        dsn <- .y
        div <- .x %>%
          # pull out needed columns
          select(marker,sample,zotu,reads) %>%
          # pivot to zotu x sample wide format
          pivot_wider(names_from="zotu",values_from="reads",values_fill = 0) %>%
          # join metadata so we can get "clean" sample names (do we need this?)
          left_join(metadata %>% select(contains("sample")), by=c("sample" = "clean_sample")) %>%
          # make sure we have the correct sample name
          select(-sample, sample = sample.y) %>%
          # order columns appropriately
          select(marker,sample,everything()) %>%
          # calculate diversity metrics
          mutate(
            shannon = diversity(pick(starts_with("taxon")),index="shannon"),
            simpson = diversity(pick(starts_with("taxon")),index="simpson"),
            richness = specnumber(pick(starts_with("taxon"))),
            trueshannon = exp(shannon),
            truesimpson = 1/(1-simpson)
          ) %>%
          # get rid of reads columns
          select(-starts_with("taxon"))
          
        # do plotting for each diversity metric
        # make a map so they can be nicely displayed
        metrics <- c(
          "Shannon" = "shannon",
          "True Shannon" = "trueshannon",
          "Simpson" = "simpson",
          "True Simpson" = "truesimpson",
          "zOTU Richness" = "richness"
        )
        # make plots for each diversity metric
        metrics %>%
          imap(~{
            # set aesthetic to div metric x marker
            ggplot(div,aes(y=marker,x=.data[[.x]])) +
              stat_summary(aes(color=marker),fun.data=data_summary,linewidth=2) +
              # geom_boxplot(aes(fill=marker)) +
              geom_point(aes(color = marker, fill = marker), shape = 21, size = 4, color = "black") +
              scale_fill_manual(values=fill_alpha(marker_pal,0.7)) + 
              xlab("Marker") + 
              ylab(.y) +
              theme_bw() +
              theme(legend.position = "none",
                    panel.grid = element_blank(),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    axis.text=element_text(size=14)) +
              coord_cartesian(xlim=c(0,max(div %>% pull(.x)))) + 
              labs(title=title_map[dsn])
          }) %>%
          # clean up the names of the returned list
          clean_names()
      })
  })

# plot true shannon across the three datasets for the unrarefied datasets
# div_plotz$raw$sharkpen$true_shannon +
  # div_plotz$raw$aquarium$true_shannon +
  # div_plotz$raw$mock$true_shannon +
  # plot_layout(axis_titles = "collect") 


# save them
if (save_pdf) {
  div_plotz %>%
    iwalk(~{
      r <- .y
      .x %>% iwalk(~{
        p <- .y
        .x %>%
          iwalk(~{
            d <- .y
            ggsave(path(fig_dir,str_glue("diversity_{p}_{d}_{r}.pdf")),.x,device=cairo_pdf,width=14,height=12,units="in")
          })
      })
    })
}

# zotu intersections (upset plots) ----------------------------------------

# intersections
# which column to show intersections for
# for reasons, these have to be expressions rather than strings
upset_cols <- c(expr(zotu),expr(family))
# label map for pretty display
ll <- c(zotu = "zOTUs", family = "families") 
# title map
pl <- c(zotu = "Intersections by zOTU", family = "Intersections by family")

upset_plotz <- datasets %>%
  map(~{
    .x %>%
      map(~{
        ds <- .x
        upset_cols %>%
          set_names() %>%
          imap(~{
            ds %>%
              # create a presence/absence dataset where we pivot to wider by marker x zotu 
              # and summarise each occurrence as a zero or one (1 = sum(reads) > 0)
              pivot_wider(id_cols={{.x}},names_from="marker",values_from="reads",values_fill = 0,values_fn=~as.integer(sum(.x) > 0)) %>%
              # plot the upset plot
              upset_plot(
                name_column={{.x}},
                data_columns=any_of(markers),
                label_top_bars = TRUE,
                label_side_bars = TRUE,
                bar_lab = str_glue("Shared {ll[as.character(.x)]}"),
                sidebar_lab = ll[as.character(.x)], 
                group_palette = marker_pal
              ) %>%
              wrap_elements() +
              labs(title=pl[as.character(.x)])
          })
      })
  })

# show upset plots for family and zotu in the unrarefied mock community
# upset_plotz$raw$mock$zotu / upset_plotz$raw$mock$family

# save them

if (save_pdf) {
  upset_plotz %>%
    iwalk(~{
      r <- .y
      .x %>%
        iwalk(~{
          p <- .y
          .x %>%
            iwalk(~{
              pl <- .y
              ggsave(path(fig_dir,str_glue("upset_{p}_{pl}_{r}.pdf")),.x,device=cairo_pdf,width=14,height=12,units="in")
            })
        })
    })
}



# expected vs unexpected species detections -------------------------------

# get ncbi ranked lineage dump
# filter it down to just bony & cartilagenous fishes for bit faster joining down below
lineage <- get_lineage(path(data_dir,"rankedlineage.dmp")) %>%
  filter(class %in% c("Actinopteri","Chondrichthyes"))


# load known invasive/introduced species
ais <- read_csv(path(data_dir,"ais.csv"))

# hacky little config map to say whether we want to see ALL possible expected species in the figure (or not)
# this is done so that the sharkpen figure doesn't show ~1200 speces in the 'expected' column
show_all <- c(
  mock = TRUE,
  aquarium = TRUE,
  sharkpen = FALSE
)

# hacky little config map to say whether we want to see invasive/introduced species
# (of course we don't care about this for the aquarium dataset)
show_inv <- c(
  mock = FALSE,
  aquarium = FALSE,
  sharkpen = TRUE
)

# whether to color invasive/introduced/whatever
color_introduced <- TRUE

# whether to do point or tile
tile <- TRUE

# create
expected_plotz <- datasets %>%
  map(~{
    .x %>%
      imap(~{
        ds <- .x
        dsn <- .y
        if (file_exists(path(data_dir,str_glue("{dsn}_species.csv")))) {
          # read expected species from file and join in lineage info from ncbi taxonomy
          # we use the lineage info to sort by family in addition to genus/species
          spp <- read_csv(path(data_dir,str_glue("{dsn}_species.csv")),col_types=cols()) %>%
            left_join(lineage,by=c("species" = "taxon")) %>%
            select(domain:genus,species) %>%
            arrange(across(domain:species))
          
          # prepare list of detected species
          sp_id <- ds %>%
            # drop unidentified things
            filter(species != "unidentified") %>%
            # group by marker and species
            group_by(marker,across(domain:species)) %>%
            summarise(reads = sum(reads)) %>%
            ungroup() %>%
            # keep only family to species
            select(-c(domain:order)) %>%
            # switch to wide (marker x detected)
            pivot_wider(names_from="marker",values_from=reads,values_fn=~.x > 0,values_fill=FALSE) 
          
          # glue text for invasive species/genera
          inv_spp <- '{species}<sup>**†**</sup>'
          inv_gen <- '{species}<sup>**‡**</sup>'
          nat_gen <- '{species}<sup>**\\***</sup>'
          # inv_spp <- '<span style="color: #C1666B">**{species}**</span>'
          # inv_gen <- '<span style="color: #D4B483">**{species}**</span>'
          # nat_gen <- '<span style="color: #4357AD">**{species}**</span>'
          # join everything together into the expected/unexpected table
          ss <- spp %>%
            # keep only family to species
            select(-c(domain:order)) %>%
            # first set all expected species to expected
            mutate(expected = TRUE) %>%
            # join in species we actually detected
            full_join(sp_id,by=c("family","genus","species")) %>%
            mutate(
              # any NA expected value is an unexpected species
              expected = replace_na(expected,FALSE),
              # any NA marker value is a negative detection
              across(any_of(markers),~replace_na(.x,FALSE))
            ) %>%
            # were there detections for any marker?
            mutate(`Any marker` = rowSums(pick(any_of(markers))) > 0) %>%
            # get rid of undetected taxa, if that's what we want to do
            filter(show_all[dsn] | `Any marker`) %>%
            # drop the any 'marker' if we're not showing everything
            { if(!show_all[dsn]) select(.,-`Any marker`) else . } %>%
            # switch to long format
            pivot_longer(-c(family:expected),names_to="marker",values_to="present")  %>%
            # drop markers with zero detections
            group_by(marker) %>%
            filter(sum(present) > 0) %>%
            ungroup() %>%
            # filter out unexpected 'Any' markers
            filter(!(!expected & marker == "Any marker")) %>%
            # create columns for plotting
            mutate( exp = if_else(expected,"Expected species","Unexpected species") ) %>%
            # smash together marker and present true/false
            unite("clr",marker,present,remove = FALSE) %>%
            { if(show_all[dsn]) mutate(.,marker = fct_relevel(marker,"Any marker",after=Inf)) else . } %>%
            # annotate known introduced/invasive species
            mutate( 
              species = case_when(
                color_introduced & 
                  show_inv[dsn] &
                  species %in% ais$scientific_name ~
                  str_glue(inv_spp),
                color_introduced & 
                  show_inv[dsn] & 
                  !expected &
                  !(species %in% ais$scientific_name) &
                  genus %in% ais$genus &
                  !(genus %in% spp$genus) ~ str_glue(inv_gen),
                color_introduced & 
                  show_inv[dsn] & 
                  !expected &
                  !(species %in% ais$scientific_name) &
                  genus %in% ais$genus &
                  genus %in% spp$genus ~ str_glue(nat_gen),
                .default = species
              )
            ) %>%
            # sort by expected and taxonomy
            arrange(desc(expected),across(family:species)) %>%
            # add row number to sort species names
            mutate(n=row_number()) %>%
            # order species names as factor levels so they display in the order we want
            mutate(species = fct_reorder(species,-n)) %>%
            # get rid of temp column
            select(-n)  
          
          # make a color palette so that each present marker gets its own
          # color, but all absent markers are just white
          mp <- c(
            c(marker_pal,"black") %>% set_names(c(str_c(names(marker_pal),"_TRUE"),"Any marker_TRUE")),
            rep("white",7) %>% set_names(str_c(c(markers,"Any marker"),"_FALSE"))
          )
          
          
          # do plot
          ggplot(ss) + 
            # points for presesence/absence of species x marker
            { if (tile) geom_tile(aes(x=marker,y=species,fill=clr),color="black") else geom_point(aes(x=marker,y=species,fill=clr),shape=21,color="black",size=2.5) } +
            # fill them by marker color
            scale_fill_manual(values=mp,guide="none",na.value = "white") + 
            # put the x axis labels on top
            scale_x_discrete(position="top") +
            # facet by expected
            facet_wrap(~exp,scales = "free") +
            # make it look nicer
            theme_bw() + 
            theme(
              strip.background = element_rect(fill="#eeeeee"),
              strip.text = element_text(face = "bold"),
              axis.title = element_blank(),
              axis.text.y = element_markdown()
            ) + 
            labs(title=title_map[dsn])
        }
      })
  })

# show a couple example plots
# expected_plotz$raw$mock
# expected_plotz$raw$aquarium
# expected_plotz$raw$sharkpen

# save expected plot figures
if (save_pdf) {
  expected_plotz %>%
    iwalk(~{
      r <- .y
      .x %>% iwalk(~{
        p <- .y
        ggsave(path(fig_dir,str_glue("expected_{p}_{r}.pdf")),.x,device=cairo_pdf,width=14,height=12,units="in")
      })
    })
}

# community multivariate statistics and pcoa plots ------------------------

community_stats <- datasets %>%
  map(~{
    dataset <- .x
    dataset %>%
      imap(~{
        dsn <- .y
        div <- .x %>%
          # pull out needed columns
          select(marker,sample,zotu,reads) %>%
          # pivot to zotu x sample wide format
          pivot_wider(names_from="zotu",values_from="reads",values_fill = 0) %>%
          # join metadata so we can get "clean" sample names (do we need this?)
          left_join(metadata %>% select(contains("sample")), by=c("sample" = "clean_sample")) %>%
          # make sure we have the correct sample name
          select(-sample, sample = sample.y) %>%
          # order columns appropriately
          select(marker,sample,everything())
        dr <- div %>%
          select(-marker) %>%
          column_to_rownames("sample") %>%
          decostand("total")
        
        sd <- div %>%
          select(sample,marker)
        
        dd <- vegdist(dr,method="bray")
        
        list(
          permanova = adonis2(dr ~ marker,data=sd,method = "bray"),
          pairwise = pairwise_adonis(dr,sd$marker,method="bray"),
          permdisp = betadisper(dd,sd$marker,bias.adjust = TRUE),
          dist = dd
        )
      })
  })

# map through permdisp objects
pcoa_plotz <- community_stats %>%
  imap(~{
    ds <- .y
    .x %>% 
      imap(~{
        dsn <- .y
        ss <- .x
        
        # get sample data
        sd <- datasets[[ds]][[dsn]] %>%
          # get distinct sample/marker combos
          distinct(sample,marker) %>%
          # make sure we have our "clean" sample names as well
          inner_join(metadata %>% select(sample,clean_sample),by=c("sample" = "clean_sample")) %>%
          # retain the correct columns
          select(-sample, sample = sample.y) %>%
          # make sure we retain missing markers so they all show up in the legends
          mutate(marker = factor(marker,levels=markers))
        
        # calculate axis explained variance
        axes <- ss$permdisp$eig / sum(ss$permdisp$eig)
        
        # make tibble of pcoa coordinates and join in sample data
        bb <- ss$permdisp$vectors %>%
          as_tibble(rownames="sample") %>%
          left_join(sd,by="sample")
        
        # plot pcoa
        ggplot(bb,aes(x=PCoA1,y=PCoA2,fill=marker)) + 
          # show.legend=TRUE makes it so missing factor levels still show up in the key
          geom_point(shape=21,size=4,color="black",show.legend=TRUE) + 
          # drop=FALSE also ensures missing factor levels are retained
          scale_fill_manual(values=marker_pal,name="Marker",drop=FALSE) +
          theme_bw() + 
          xlab(str_glue("PCOA1 ({scales::percent(axes[1],accuracy=0.1)})")) + 
          ylab(str_glue("PCOA2 ({scales::percent(axes[2],accuracy=0.1)})")) + 
          labs(title=title_map[dsn])
      })
  })

# make composite figures for each major dataset (raw vs rarefied)
pcoa_composites <- pcoa_plotz %>%
  map(~{
    .x %>%
      reduce(`+`)
  })

# show plots, smashing together similar legends
# pcoa_composites$raw + plot_layout(guides="collect")
# pcoa_composites$rarefied + plot_layout(guides="collect")

# save them
if (save_pdf) {
  pcoa_plotz %>%
    iwalk(~{
      r <- .y
      .x %>%
        iwalk(~{
          p <- .y
          ggsave(path(fig_dir,str_glue("pcoa_{p}_{r}.pdf")),.x,device=cairo_pdf,width=14,height=12,units="in")
        })
    })
}



# cluster plots -----------------------------------------------------------

# use previously-calculated distance matrices do do cluster plots
# this'll spit some warnings because the new ggplot.ggdendro uses some
# old-style code
cluster_plotz <- community_stats %>%
  imap(~{
    ds <- .y
    .x %>%
      imap(~{
        dsn <- .y
        dd <- .x$dist
        
        # do a cluster analysis and convert to ggdend object
        upgma <- dd %>% 
          hclust(method = "average") %>%
          as.dendrogram() %>%
          as.ggdend()
        
        # plot it using ggplot.ggdend
        uplot <- ggplot(upgma,hang=-1,angle=45) +
          expand_limits(y=-0.15,x=-0.3)
        
        # do a cluster analysis and convert to ggdend object
        ward <- dd %>%
          hclust(method = "ward.D2") %>%
          as.dendrogram() %>%
          as.ggdend()
        
        # plot it using ggplot.ggdend
        wplot <- ggplot(ward,hang=-1,angle=45) +
          expand_limits(y=-0.3,x=-0.3)
        
        # return list of both plots
        list(upgma = uplot, ward = wplot)
        
      })
  })

# make composites of the ward clusters
ward_composite <- cluster_plotz %>% 
  map(~{
    .x %>% 
      imap(~.x$ward + labs(title=title_map[.y])) %>%
      reduce(`+`)
  })
# ward_composite$raw


# make composites of the upgma clusters
upgma_composite <- cluster_plotz %>% 
  map(~{
    .x %>% 
      imap(~.x$upgma + labs(title=title_map[.y])) %>%
      reduce(`+`)
  })
# upgma_composite$raw


# save them
if (save_pdf) {
  cluster_plotz %>%
    iwalk(~{
      r <- .y
      .x %>%
        iwalk(~{
          p <- .y
          .x %>%
            iwalk(~{
              pl <- .y
              ggsave(path(fig_dir,str_glue("cluster_{p}_{pl}_{r}.pdf")),.x,device=cairo_pdf,width=14,height=12,units="in")
            })
        })
    })
}
