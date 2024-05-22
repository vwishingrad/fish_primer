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
library(circlize)
library(ggplotify)


# setup -------------------------------------------------------------------

# whether to save to pdf all the time
save_pdf <- FALSE

# prefill list of plot tags
plot_tags <- list(str_glue("({letters})"))

# source utility functions
source("util.R")

# make figures directory
fig_dir <- here("output","figures")
dir_create(fig_dir,recurse = TRUE)

# make tables directory
tbl_dir <- here("output","tables")
dir_create(tbl_dir,recurse = TRUE)

# get rid of annoying '`summarise()` has grouped output by' message
options(dplyr.summarise.inform = FALSE)

# get base data directory
data_dir <- here("data")

# read markers
markers <- read_lines(path(data_dir,"markers.txt"))

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

# whether to rarefy at all
rarefy <- TRUE

# number of rarefaction permutations
rarefy_perm <- 100

# whether to filter out unidentified taxa
filter_unid <- TRUE

# helper to relevel taxonomic factors
relevel_unid <- function(f) {
  unid <- grep("unidentified",f,value=TRUE) %>%
    unique() %>%
    sort()
  fct_relevel(f,unid,after=Inf)
}

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
fishes_raw <- markers %>%
  map(~{
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
        zotus = str_glue("{marker}({str_c(otu,collapse=',')})"), 
        zotu_count = n(),
        marker = marker
      ) %>%
      ungroup() %>%
      select(domain:species,marker,representative,zotus,zotu_count,all_of(samples))
  }) %>%
  list_rbind() %>%
  # these operations are done on the fully concatenated dataset
  mutate(
    # replace NA taxa with blank string
    across(domain:species,~replace_na(.x,"")),
    # replace NA reads with zeroes
    across(where(is.numeric),~replace_na(.x,0))
  )

fishes <- fishes_raw %>%
  
  # this is a convoluted way to get zeroes in the per-marker zotu count
  select(domain:species,marker,zotu_count) %>%
  pivot_wider(names_from = "marker",values_from="zotu_count",values_fill=0) %>%
  pivot_longer(all_of(markers),names_to = "marker",values_to = "zotu_count") %>%
  full_join(fishes_raw,by=c("domain","kingdom","phylum","class","order","family","genus","species","marker")) %>%
  mutate(zotu_count = coalesce(zotu_count.x,zotu_count.y)) %>%
  select(all_of(names(fishes_raw))) %>%
  arrange(across(domain:species),marker) %>%
  mutate(across(all_of(all_samples),~replace_na(.x,0))) %>%
  replace_na(list(representative="",zotus="")) %>%

  # sum reads by taxon and concatenate zotus/representatives
  group_by(across(domain:species)) %>%
  summarise(
    # concatenate representative zotus
    representative=str_c(representative,collapse=','),
    # concatenate zotu lists
    zotus=str_c(zotus,collapse=','),
    total_zotu_count = sum(zotu_count),
    marker_zotu_count = str_c(str_glue("{marker}:{zotu_count}"),collapse=","),
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
  # get rid of samples with total read counts under minimum threshold
  select(-all_of(samples),where(~is.numeric(.x) && sum(.x) >= min_total)) %>%
  
  # here we try to get clever about how to deal with blank taxa
  # anything that's a blank between non-blanks (e.g., the missing order of pomacentrids)
  # gets labeled as "incertae sedis". anything that's blank at the end gets labeled
  # as "unidentified" and then the last identified taxon
  # first, pivot to long format
  pivot_longer(domain:species,names_to="taxon_level",values_to="taxon") %>%
  select(taxon_level,taxon,everything()) %>%
  # now group by individual taxon
  group_by(zotu) %>%
  group_modify(~{
    # are any levels blank?
    if (any(.x$taxon == "")) {
      # which ones are blank?
      blanks <- which(.x$taxon == "") 
      # are the places before the blanks also blank?
      before_blanks <- replace_na(.x$taxon[blanks-1] == "",TRUE)
      # this is a special case where the first slot (domain) is blank
      if (any(blanks == 1)) before_blanks <- c(FALSE,before_blanks)
      # are the palces after the blanks also blank?
      after_blanks <- replace_na(.x$taxon[blanks+1] == "",TRUE)
      # anything blank that's sandwiched between non-blanks gets 'incertae sedis'
      .x$taxon[blanks[which(!before_blanks & !after_blanks)]] <- "incertae sedis"
      # do we still have blanks?
      # if we do, I think we can assume they're at the end 
      # (tho I guess two consecutive blanks would mess this up)
      if (any(.x$taxon == "")) {
        # get last non-blank slot 
        last_nonblank <- min(which(.x$taxon == ""))-1
        # get last non-blank name (if it's the first slot, just leave it blank)
        last_name <- ifelse(last_nonblank > 0,str_c(" ",.x$taxon[last_nonblank]),"")
        # make the new name "unidentified <last non-blank name>"
        unid <- str_glue("unidentified{last_name}")
        # set all the remaining blanks to that new name
        .x$taxon <- replace(.x$taxon,which(.x$taxon == ""),unid)
      }
    }
    return(.x)
  }) %>%
  # ungroup or all sorts of things get screwy
  ungroup() %>%
  # go back to wide
  pivot_wider(names_from="taxon_level",values_from="taxon") %>%
  mutate(
    # in the special case where species is "unidentified <genus>",
    # let's change it to "<genus> sp."
    species = case_when(
      species == str_glue("unidentified {genus}") ~ str_glue("{genus} sp."),
      .default = species
    )
  ) %>%
  select(domain:species,everything()) %>%
  # once more, get rid of taxa with zero reads
  mutate(
    tot = rowSums(pick(matches(unlist(dataset_map))))
  ) %>%
  filter(tot > 0) %>%
  select(-tot) %>%
  # this smashes all the "unidentified" taxa to the end of the factor-level order
  mutate( across(domain:species,~relevel_unid(.x)) ) %>%
  arrange(class,family) %>% 
  mutate(
    # order families' factor levels alphabetically by class first, then by family name
    family = fct_reorder(family,order(class)),
  )


rr <- {if (rarefy) c(FALSE,TRUE) else c(FALSE)}
rn <- {if (rarefy) c("raw","rarefied") else c("raw")}

# use the map to break out the different sample types and pivot to long format
# do one that's rarefied and one that's not
datasets <- rr %>%
  set_names(rn) %>%
  map(~{
    rarefy <- .x
    dataset_map %>%
      map(~{
        ff <- fishes_filtered %>%
          # grab samples that match the current sample type
          select(domain:species,contains("zotu"),matches(.x))
        # rarefy datasets if so desired
        if (rarefy) {
          # pull out taxonomy part to stick back on later
          taxonomy <- ff %>% 
            select(contains("zotu"),domain:species) 
          ff <- ff %>%
            # get only zotu and samples
            select(zotu,any_of(samples)) %>%
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
            select(domain:species,contains("zotu"),everything())
        }
        ff %>%
          # switch to long format
          pivot_longer(-c(domain:species,contains("zotu")),names_to="sample",values_to="reads") %>%
          # join in sample metadata
          left_join(metadata,by=c("sample" = "clean_sample"),suffix=c("","_display")) %>%
          # group by zotu (which I guess it probably is by default)
          group_by(domain,kingdom,phylum,class,order,family,genus,species,marker,type,sample,across(contains("zotu"))) %>%
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
            marker = factor(marker,levels=markers),
            marker_zotu_count = str_extract(marker_zotu_count,str_glue("(?<={marker}:)[0-9]+")) %>%
              as.numeric()
          )
      })
  })

# get counts of unidentified taxa
unid_counts <- datasets %>%
  map(~{
    .x %>%
      imap(~{
        .x %>%
          distinct(across(domain:species)) %>%
          filter(if_any(domain:species,~str_detect(.x,"unidentified") | str_detect(.x,"sp\\.$"))) %>%
          nrow()
      })
  })


all_taxa <- datasets$raw %>%
  map(~.x %>% select(domain:species)) %>%
  list_rbind() %>%
  # unlist() %>%
  distinct(domain,kingdom,phylum,class,order,family,genus,species,.keep_all = TRUE) %>%
  arrange(class,family,species) 

# make a consistent color palette for markers
marker_pal <- paletteer_d("ggthemes::Tableau_10",n=length(markers)) %>%
  as.character() %>%
  set_names(markers)

# Let's make standardized color palettes for the various groups
palettes <- c("class","order","family","species") %>%
  set_names() %>%
  map(~{
    feesh <- all_taxa %>%
      filter(class == "Actinopteri") %>%
      distinct(.data[[.x]],.keep_all = TRUE) %>%
      pull(.data[[.x]])    
    sharks <- all_taxa %>%
      filter(class == "Chondrichthyes") %>%
      distinct(.data[[.x]],.keep_all = TRUE) %>%
      pull(.data[[.x]])    
    
    c(
      colorRampPalette(brewer.pal(12,"Paired"))(length(feesh)),
      colorRampPalette(brewer.pal(8,"Set2"))(length(sharks))
    ) %>% 
      set_names(c(feesh,sharks))
  })


# start here


# generate tables ---------------------------------------------------------
datasets$raw %>%
  iwalk(~{
    .x <- .x %>%
      filter(if_all(domain:genus,~!str_detect(.x,"unidentified")) & !str_detect(species,"sp\\.$")) %>%
      rename(Marker=marker)
    ds <- .x %>%
      mutate(Marker = fct_expand(Marker,"Overall")) %>%
      group_by(Marker) %>%
      summarise(
        `Unique families` = n_distinct(family),
        `Unique species` = n_distinct(zotu),
        `zOTUs` = sum(marker_zotu_count[!duplicated(zotu)])
      ) %>%
      ungroup() %>%
      rbind(
        list(
          Marker = "Overall",
          `Unique families` = n_distinct(.x$family),
          `Unique species` = n_distinct(.x$zotu),
          `zOTUs` = sum(.x$marker_zotu_count[!duplicated(.x$zotu)]) 
        )
      )
    write_csv(ds,path(tbl_dir,str_glue("{.y}_summary.csv")))
  })

datasets$raw %>%
  imap(~{
    .x %>%
      rename(Marker = marker) %>%
      mutate(Marker = fct_expand(Marker,"Overall")) %>%
      group_by(Marker) %>%
      summarise(
        `Unique families` = n_distinct(family),
        `Unique species` = n_distinct(zotu),
        `zOTUs` = sum(marker_zotu_count[!duplicated(zotu)])
      ) %>%
      ungroup() %>%
      rbind(
        list(
          Marker = "Overall",
          `Unique families` = n_distinct(.x$family),
          `Unique species` = n_distinct(.x$zotu),
          `zOTUs` = sum(.x$marker_zotu_count[!duplicated(.x$zotu)]) 
        )
      ) %>%
      mutate(`Sample type` = title_map[.y]) 
  }) %>%
  list_rbind() %>%
  write_csv(path(tbl_dir,"all_summary.csv"))


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
        # map through datasets
        dataset %>%
          imap(~{
            # save dataset name
            nn <- .y
            # sum up relative reads for marker/plot level combos
            dd <- .x %>%
              filter(!filter_unid | .data[[pl]] != "unidentified") %>%
              { if (filter_unid) mutate(.,"{pl}" := fct_recode(.data[[pl]],NULL="unidentified")) else . } %>%
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
              geom_col(aes(x=marker,y=rel,fill=.data[[pl]]),show.legend = TRUE) +
              scale_x_discrete(drop = FALSE) +
              scale_y_continuous(expand = expansion(mult=c(0.01,NA)) ) +
              scale_fill_manual(
                values=palettes[[pl]],
                drop=FALSE,
                name=str_to_sentence(pl)
              ) + 
              theme_bw() + 
              theme(
                panel.border = element_blank(),
                axis.line = element_line(color="black"),
                panel.grid = element_blank(),
                axis.text.x = element_text(angle=90),
                axis.title = element_text(face="bold")
              ) + 
              labs(y=str_glue("Relative sequence abundance"),x="Marker")
          })
      }) 
  })

# rel_taxon_plotz$raw$family$aquarium
# rel_taxon_plotz$raw$family$mock
# rel_taxon_plotz$raw$family$sharkpen

# reduce all these plots to one big one using patchwork

# make composite plots for each plot level
rel_taxon_composites <- rel_taxon_plotz %>%
  map(~{
    .x %>%
      imap(~{
        .x %>%
          map(~.x + theme(axis.text.x = element_text(angle=45,hjust=1))) %>%
          reduce(`/`) +
          plot_annotation(tag_levels = plot_tags) + 
          plot_layout(guides="collect", axis_titles = "collect")
      })
  })

if (save_pdf) {
  rel_taxon_composites %>%
    iwalk(~{
      r <- .y
      .x %>%
        iwalk(~{
          pl <- .y
          .x <- .x &
            guides(fill = guide_legend(ncol = 3)) & 
            theme(
              axis.text = element_text(size=20),
              legend.text = element_text(size=24),
              axis.title = element_text(size=24),
              axis.title.y = element_text(size=24),
              legend.title = element_text(size=24),
              legend.key.size = unit(26,units="pt"),
              plot.tag = element_text(size=30)
            )
          ggsave(path(fig_dir,str_glue("taxon_composite_{pl}_{r}.pdf")),.x,device=cairo_pdf,width=24,height=28,units="in")
        })
    })
}
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
        dataset %>%
          imap(~{
            nn <- .y
            # calculate relative zotu abundance (not read abundance)
            # for marker and plot level
            dd <- .x %>%
              filter(reads > 0) %>%
              filter(!filter_unid | .data[[pl]] != "unidentified") %>%
              group_by(marker,across(all_of(pl))) %>%
              summarise(n = n_distinct(zotu)) %>%
              group_by(marker) %>%
              mutate(rel = n / sum(n)) %>%
              ungroup()
            if (!zotu_rel)
              dd <- dd %>% mutate(rel=n)
            # do plotting
            ggplot(dd) + 
              geom_col(aes(x=marker,y=rel,fill=.data[[pl]]),show.legend = TRUE) + 
              scale_x_discrete(drop = FALSE) +
              # scale_y_continuous(expand = c(0, 0), limits = c(0,1.01) )+
              scale_fill_manual(values=palettes[[pl]],name=str_to_sentence(pl),drop=FALSE) + 
              scale_y_continuous(expand = expansion(mult=c(0.01,NA)) ) +
              theme_bw() + 
              theme(
                panel.border = element_blank(),
                axis.line = element_line(color="black"),
                panel.grid = element_blank(),
                axis.text.x = element_text(angle=90),
                axis.title = element_text(face="bold")
              ) + 
              labs(y="Number of unique taxa",x="Marker")
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
          reduce(`/`) +
          plot_layout(axis_titles ="collect",guides="collect") + 
          plot_annotation(tag_levels = plot_tags)
      })
  })

if (save_pdf) {
  zotu_composites %>%
    iwalk(~{
      r <- .y
      .x %>%
        iwalk(~{
          pl <- .y
          .x <- .x &
            guides(fill = guide_legend(ncol = 3)) & 
            theme(
              axis.text = element_text(size=20),
              legend.text = element_text(size=24),
              axis.title = element_text(size=24),
              axis.title.y = element_text(size=24),
              legend.title = element_text(size=24),
              legend.key.size = unit(26,units="pt"),
              plot.tag = element_text(size=30)
            )
          ggsave(path(fig_dir,str_glue("zotu_composite_{pl}_{r}.pdf")),.x,device=cairo_pdf,width=24,height=28,units="in")
        })
    })
}
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
              geom_tile(aes(x=marker,y=.data[[pl]],fill=rel),color="grey8") + 
              scale_fill_paletteer_c("viridis::turbo",name="Relative\nabundance")  +
              scale_y_discrete(limits=rev) + 
              labs(x="Marker",y=str_to_sentence(pl)) + 
              theme(axis.title = element_text(face="bold"))
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

# make composite plots for each
zotu_heats <- rel_zotu_plotz %>%
  map(~{
    .x %>%
      imap(~{
        .x %>%
          map(~.x + theme(axis.text.x = element_text(angle=45,hjust=1))) %>%
          reduce(`+`) +
          plot_layout(axis_titles = "collect") + 
          plot_annotation(title="",tag_levels = plot_tags) &
          theme(plot.tag = element_text(size=18))
      })
  })

if (save_pdf) {
  zotu_heats %>%
    iwalk(~{
      r <- .y
      .x %>%
        iwalk(~{
          pl <- .y
          ggsave(path(fig_dir,str_glue("zotu_heatmap_{pl}_{r}.pdf")),.x,device=cairo_pdf,width=18,height=8,units="in")
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
upset_cols <- c(expr(family),expr(zotu))
# label map for pretty display
ll <- c(zotu = "Taxa", family = "Families", species = "Species") 
# title map
pl <- c(zotu = "Intersections by taxon", family = "Intersections by family", species = "Intersections by species")

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
              wrap_elements() 
          })
      })
  })

# show upset plots for family and zotu in the unrarefied mock community
# upset_plotz$raw$sharkpen$family / upset_plotz$raw$sharkpen$zotu + plot_annotation(tag_levels = plot_tags)
# upset_plotz$raw$aquarium$family / upset_plotz$raw$aquarium$zotu + plot_annotation(tag_levels = plot_tags)
# upset_plotz$raw$mock$family / upset_plotz$raw$mock$zotu + plot_annotation(tag_levels = plot_tags)

# put the plots together for each sample type
upset_composites <- upset_plotz %>%
  map(~{
    .x %>%
      map(~{
        .x %>%
          reduce(`/`) + 
          plot_annotation(tag_levels = plot_tags)
      })
  })

# save them
if (save_pdf) {
  upset_composites %>%
    iwalk(~{
      r <- .y
      .x %>%
        iwalk(~{
          p <- .y
          ggsave(path(fig_dir,str_glue("upset_{p}_{r}.pdf")),.x,device=cairo_pdf,width=15,height=11,units="in")
        })
    })
}



# expected vs unexpected species detections -------------------------------

# get ncbi ranked lineage dump
# filter it down to just bony & cartilagenous fishes for bit faster joining down below
# actually, forget about it. we just put the info straight into the expect species spreadsheets
# lineage <- get_lineage(path(data_dir,"rankedlineage.dmp")) %>%
#   filter(class %in% c("Actinopteri","Chondrichthyes"))


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

exp_map = list(
  sharkpen = c("Not known from Hawai‘i","Known from Hawai‘i"),
  aquarium = c("Not known from aquarium tank","In aquarium tank"),
  mock = c("Not in mock community","In mock community")
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
          spp <- read_csv(path(data_dir,str_glue("{dsn}_species.csv")),col_types=cols())
          
          # prepare list of detected species
          sp_id <- ds %>%
            # drop unidentified things
            # filter(species != "unidentified") %>%
            filter(
              !str_detect(species,"sp\\.$"),
              !str_detect(species,"^unidentified")
            ) %>%
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
            mutate( exp = exp_map[[dsn]][as.integer(expected)+1] ) %>%
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
              ),
              species = str_glue("*{species}*")
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
            # points/tiles for presesence/absence of species x marker
            { if (tile) geom_tile(aes(x=marker,y=species,fill=clr),color="gray10") else geom_point(aes(x=marker,y=species,fill=clr),shape=21,color="black",size=2.5) } +
            # fill them by marker color
            scale_fill_manual(values=mp,guide="none",na.value = "white") + 
            # put the x axis labels on top
            scale_x_discrete(position="top") +
            # facet by expected
            ggforce::facet_col(~exp,scales="free",space="free") + 
            # make it look nicer
            theme_bw() + 
            theme(
              strip.background = element_rect(fill="#eeeeee"),
              strip.text = element_text(face = "bold"),
              axis.title = element_blank(),
              axis.text.y = element_markdown(),
              axis.text.x = element_text(face="bold")
            ) 
        }
      })
  })

# show a couple example plots
# expected_plotz$raw$mock
# expected_plotz$raw$aquarium
# expected_plotz$raw$sharkpen

# another hacky little map so we get the right plot sizes
plotsizes <- list(
  sharkpen = c(7.7,7.7),
  aquarium = c(7.7,7.7),
  mock = c(9.8,10.7)
)

# save expected plot figures
if (save_pdf) {
  expected_plotz %>%
    iwalk(~{
      r <- .y
      .x %>% iwalk(~{
        p <- .y
        ggsave(path(fig_dir,str_glue("expected_{p}_{r}.pdf")),.x,device=cairo_pdf,width=plotsizes[[p]][1],height=plotsizes[[p]][2],units="in")
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
        
        dd <- vegdist(dr,method="jaccard")
        
        list(
          permanova = adonis2(dr ~ marker,data=sd,method = "jaccard"),
          pairwise = pairwise_adonis(dr,sd$marker,method="jaccard"),
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
          geom_point(shape=21,size=4,color="black",show.legend=TRUE) + 
          scale_fill_manual(values=marker_pal,name="Marker",drop=FALSE) +
          theme_bw() + 
          xlab(str_glue("PCoA1 ({scales::percent(axes[1],accuracy=0.1)})")) + 
          ylab(str_glue("PCoA2 ({scales::percent(axes[2],accuracy=0.1)})"))
      })
  })

# make composite figures for each major dataset (raw vs rarefied)
pcoa_composites <- pcoa_plotz %>%
  map(~{
    .x %>%
      reduce(`+`) + 
      plot_layout(axes = "collect", guides = "collect") + 
      plot_annotation(tag_levels = plot_tags) &
      theme(plot.tag = element_text(face="bold"))
  })


if (save_pdf) {
  pcoa_composites %>%
    iwalk(~{
      r <- .y
      ggsave(path(fig_dir,str_glue("pcoa_composite_{r}.pdf")),.x,device=cairo_pdf,width=12,height=4,units="in")
    })
}

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


# chord plots -------------------------------------------------------------
make_chord <- function(dd,margin,pal,units="",group=NULL) {
  circos.clear()
  if (units == "in") {
    margin <- inches_h(margin)
  }
  circos.par(circle.margin = margin)
  chordDiagram(dd, annotationTrack = "grid", preAllocateTracks = 1,grid.col = pal,group=group)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    if (CELL_META$sector.index %in% markers) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1]+CELL_META$cell.height, CELL_META$sector.index,
                  facing = "bending.inside", adj = c(0.5,1),cex=1.2)
    } else {
      cex = if_else(CELL_META$cell.width < 2,0.9,1)
      yoffs = 0
      circos.text(CELL_META$xcenter, CELL_META$ylim[1]+yoffs, CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=cex)
    }
  }, bg.border = NA)
  circos.clear()
}


chord_plotz <- datasets$raw %>%
  imap(~{
    dsn <- .y
    dd <- .x %>%
      filter(family != "unidentified") %>%
      count(marker,class,family) %>%
      mutate(
        marker = factor(marker,levels=sort(markers)),
        marker = fct_relevel(marker,"MiFish_U","MiFish_E",after=Inf)
      ) %>%
      arrange(class,family,marker) %>%
      select(marker,family,n) 
    
    
    mk <- dd %>%
      distinct(marker) %>%
      pull(marker) %>%
      as.character()
    mk <- rep("marker",length(mk)) %>%
      set_names(mk)
    fm <- .x %>%
      arrange(class,desc(family)) %>%
      distinct(family,.keep_all = TRUE) %>%
      filter(family %in% dd$family) %>%
      select(family,class) %>%
      mutate(class = if_else(class == "Chondrichthyes","sharks","fishes"),family=as.character(family)) %>%
      deframe()
    grp <- c(mk,fm)
    as.ggplot(function() make_chord(dd,c(0.000000001,1,0.000000001,0.000000001),c(marker_pal,palettes$family),units="in",group=grp))
  })

# chord_plotz$mock

if (save_pdf) {
  chord_plotz %>%
    iwalk(~{
      ggsave(path(fig_dir,str_glue("chord_{.y}.pdf")),.x,device=cairo_pdf,width=10.5,height=10.5,units="in")
    })
}

# raw data supplemental bar plots -----------------------------------------

raw_plot_data <- markers %>%
  set_names() %>%
  map(~{
    # save marker name into something nicer than `.x`
    marker <- .x
    
    # read marker data and filter to just chordates
    ds <- read_tsv(path(data_dir,str_glue("{marker}_data.tsv")),col_types = cols()) %>%
      mutate(
        # replace "dropped" with blank string
        across(class:species,~replace(.x,which(.x == "dropped"),"")),
        # sum across blanks into column `blanks`
        blanks = rowSums(pick(matches("Blank")),na.rm=TRUE),
        # subtract blank reads from sample reads
        across(-c(domain:seq_length,blanks),~.x-blanks),
        # set negative reads to zero
        across(where(is.numeric),~replace(.x,which(.x < 0),0))
      ) %>%
      select(-starts_with(marker),where(~is.numeric(.x) && sum(.x) > 0),-blanks) %>%
      pivot_longer(-c(domain:seq_length),names_to = "sample",values_to="reads") %>%
      mutate(
        grouping = case_when(
          domain == "Bacteria" ~ "Bacteria",
          class == "Mammalia" ~ "Mammals",
          class == "Chondrichthyes" ~ "Sharks & rays",
          class == "Actinopteri" ~ "Bony fishes",
          .default = "Other"
        ),
        grouping = factor(grouping,levels=c("Bacteria","Other","Mammals","Sharks & rays","Bony fishes"))
      ) %>%
      filter(reads > 0)  %>%
      mutate(
        sample_type = case_when(
          str_detect(sample,regex(dataset_map$sharkpen,ignore_case = TRUE)) ~ "Shark pen",
          str_detect(sample,regex(dataset_map$aquarium,ignore_case = TRUE)) ~ "Waikīkī Aquarium",
          str_detect(sample,regex(dataset_map$mock,ignore_case = TRUE)) ~ "Mock community"
        )
      ) %>%
      # group_by(sample_type,grouping) %>%
      # mutate(total_zotus = n_distinct(OTU)) %>%
      # ungroup() %>%
      group_by(sample_type,sample,grouping) %>%
      # summarise(reads = sum(reads), zotus = n_distinct(OTU), total_zotus = unique(total_zotus)) %>%
      summarise(reads = sum(reads), zotus = n_distinct(OTU)) %>%
      ungroup() %>%
      group_by(sample) %>%
      mutate(rel_reads = reads / sum(reads), rel_zotus = zotus/sum(zotus))  %>%
      arrange(sample_type,grouping)
  })

raw_pal <- c(
  "Bony fishes" = "#3B01FE",
  "Sharks & rays" = "#A9A9A9",
  "Mammals" = "#EE82EF",
  "Bacteria" = "#ED0105",
  "Other" = "#F2A503"
)

raw_plotz <- raw_plot_data %>%
  map(~{
    read_abundance <- ggplot(.x) + 
      geom_col(aes(x=sample,y=reads,fill=grouping)) + 
      scale_fill_manual(values = raw_pal, name = "Group") +
      scale_y_continuous(labels=scales::comma, expand = expansion(mult=c(0.01,NA))) +
      # facet_wrap(~sample_type,scales="free_x",strip.position = "bottom") +
      labs(x = "",y="Abundance")
      
    
    read_rel <- ggplot(.x) + 
      geom_col(aes(x=sample,y=rel_reads,fill=grouping)) + 
      scale_fill_manual(values = raw_pal, name = "Group") + 
      scale_y_continuous(labels=scales::percent, expand = expansion(mult=c(0.01,NA))) +
      # facet_wrap(~sample_type,scales="free_x",strip.position = "bottom") +
      labs(x = "Sequence reads",y="Relative abundance") + 
      theme(axis.title.x = element_blank())
    
    zotu_abundance <- ggplot(.x) + 
      geom_col(aes(x=sample,y=zotus,fill=grouping)) + 
      scale_fill_manual(values = raw_pal, name = "Group") + 
      scale_y_continuous(labels=scales::comma, expand = expansion(mult=c(0.01,NA))) +
      # facet_wrap(~sample_type,scales="free_x",strip.position = "bottom") +
      labs(x="",y="Abundance") + 
      theme(axis.title.x = element_blank()) 
    
    zotu_rel <- ggplot(.x) + 
      geom_col(aes(x=sample,y=rel_zotus,fill=grouping)) + 
      scale_fill_manual(values = raw_pal, name = "Group") + 
      scale_y_continuous(labels=scales::percent, expand = expansion(mult=c(0.01,NA))) +
      # facet_wrap(~sample_type,scales="free_x",strip.position = "bottom") +
      labs(x="zOTUs",y="Relative abundance") 
    
    layout <- "
      12
      34
    "
    
    read_abundance + zotu_abundance + read_rel + zotu_rel +
      plot_layout(ncol = 2, guides = "collect", axis_titles = "collect") +
      plot_annotation(tag_levels = plot_tags) &
      facet_wrap(~sample_type,scales="free_x",strip.position = "bottom") &
      theme_bw() & 
      theme(
        panel.border = element_blank(),
        axis.line = element_line(color="black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.title = element_text(face="bold"),
        strip.background = element_rect(fill=NA,color=NA)
      )
  })

if (save_pdf) {
  raw_plotz %>%
    iwalk(~{
      ggsave(path(fig_dir,str_glue("reads_zotus_{.y}.pdf")),.x,device=cairo_pdf,width=15,height=10,units="in")
    })
}

num <- function(x,a=NULL) scales::label_comma(accuracy=a)(x)

marker_summary_tables <- raw_plot_data %>%
  map(~{
    .x %>%
      mutate(grp = fct_other(grouping,drop=c("Sharks & rays","Bony fishes"),other_level = "fishes")) %>%
      group_by(sample_type,sample) %>%
      summarise(
        fish = sum(reads[grp == "fishes"]),
        reads = sum(reads),
        fish_zotus = sum(zotus[grp == "fishes"]),
        zotus = sum(zotus),
      ) %>%
      ungroup() %>%
      group_by(sample_type) %>%
      summarise(
        `Total reads` = num(sum(reads)),
        `Mean reads per sample` = str_glue("{num(mean(reads))} ± {num(sd(reads))}"),
        
        `Total fish reads` = num(sum(fish)),
        `Mean fish reads per sample` = str_glue("{num(mean(fish))} ± {num(sd(fish))}"),
        
        `Total zOTUs` = num( sum(zotus)  ),
        `Mean zOTUs per sample` = str_glue("{num(mean(zotus))} ± {num(sd(zotus))}"),
        
        `Total fish zOTUs` = num( sum(fish_zotus) ),
        `Mean fish zOTUs per sample` = str_glue("{num(mean(fish_zotus))} ± {num(sd(fish_zotus))}")
      ) %>%
      ungroup() %>%
      rename(`Sample type` = sample_type)
  })

View(marker_summary_tables$Berry_16S)

marker_summary_tables %>%
  iwalk(~{
    write_tsv(.x,path(tbl_dir,str_glue("{.y}_marker_summary.tsv")))  
  })
