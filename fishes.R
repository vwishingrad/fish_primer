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


# setup -------------------------------------------------------------------

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
get_lineage <- function() {
  nlf <- path(data_dir,"rankedlineage.dmp")
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
          mutate(across(domain:species,~replace(.x,which(.x == ""),"unidentified")))
      })
  })

# start here

# taxon bar plots ---------------------------------------------------------

# whether to filter out unidentified taxa
filter_unid <- TRUE

# taxonomic levels at which to group bar plots
plot_levels <- c("family","order","class")

# create a list of bar plots by plot level and dataset
#map through plot levels
rel_taxon_plotz <- plot_levels %>%
  set_names() %>%
  map(~{
    # save plot level into `pl`
    pl <- .x
    
    # make consistent taxon palette across datasets
    # pull all possible values for the plot level
    groups <- datasets %>%
      map(~.x %>% pull(all_of(pl))) %>%
      unlist() %>%
      sort() %>%
      unique()    
    
    # make palette generator and generate a palette
    pal = colorRampPalette(brewer.pal(12, "Paired"))(length(groups))
    # name the palette
    names(pal) <- groups
    
    # map through datasets
    datasets %>%
      imap(~{
        # save dataset name
        nn <- .y
        
        # sum up relative reads for marker/plot level combos
        dd <- .x %>%
          filter(!filter_unid | .data[[pl]] != "unidentified") %>%
          group_by(marker,across(all_of(pl))) %>%
          summarise(rel = sum(rel)) %>%
          ungroup() 
        
        # do plotting
        ggplot(dd) + 
          geom_col(aes(x=marker,y=rel,fill=.data[[pl]])) +
          scale_x_discrete(drop = FALSE) +
          scale_y_continuous(expand = c(0, 0), limits = c(0,1) )+
          scale_fill_manual(values=pal,drop=FALSE) + 
          theme_bw() +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.text=element_text(size=12),
                panel.border = element_blank(),
                axis.title.y = element_text(size=14),
                axis.line.x = element_line(color="black", linewidth = .5),
                axis.line.y = element_line(color="black", linewidth = .5)) +
          labs(title=nn)
      })
  })

# reduce all these plots to one big one using patchwork
# (it looks like shit because there's too much going on)
rel_taxon_composite <- rel_taxon_plotz %>%
  imap(~{
    .x %>%
      reduce(`+`) + 
      plot_annotation(title=.y) +
      plot_layout(guides="collect")
  }) %>%
  reduce(`/`)
rel_taxon_composite


# zotu bar plots ----------------------------------------------------------

# plot zotus by taxonomic level
plot_levels <- c("family","species")


# are zotu plots relative?
zotu_rel <- FALSE

# create a list of bar plots by plot level and dataset
# map through plot levels
zotu_plotz <- plot_levels %>%
  set_names() %>%
  map(~{
    # save plot level
    pl <- .x
    # map through datasets
    
    # make consistent taxon palette across datasets
    # pull all possible values for the plot level
    groups <- datasets %>%
      map(~.x %>% pull(all_of(pl))) %>%
      unlist() %>%
      sort() %>%
      unique()    
        # make palette generator and generate a palette
    pal = colorRampPalette(brewer.pal(12, "Paired"))(length(groups))
    # name the palette
    names(pal) <- groups
    
    datasets %>%
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
          scale_fill_manual(values=pal) + 
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
          labs(title=nn)
      })
  })

# we're not gonna make a composite with these because the legends make them out of control
# but here are examples excluding the legends for family:
zotu_plotz$family$sharkpen + zotu_plotz$family$aquarium + zotu_plotz$family$mock & theme(legend.position = "none")
# and species:
zotu_plotz$species$sharkpen + zotu_plotz$species$aquarium + zotu_plotz$species$mock & theme(legend.position = "none")

# diversity plots ---------------------------------------------------------

# create box plots for diversity metrics

# map through datasets in the original wide format
div_plotz <- dataset_map %>%
  imap(~{
    # dataset name
    dsn <- .y
    div <- fishes_filtered %>%
      # match current dataset
      select(zotu,matches(.x)) %>%
      # make it longer
      pivot_longer(-zotu,names_to="sample",values_to="reads") %>%
      # join in metadata
      left_join(metadata,by=c("sample" = "clean_sample"),suffix=c("","_display")) %>%
      # grab columns we care about
      select(zotu,marker,sample=sample_display,reads) %>%
      # make it wider with zotus as the columns, fill NAs with zeroes
      pivot_wider(names_from="zotu",values_from="reads",values_fill = 0) %>%
      # calculate diversity metrics
      mutate(
        shannon = diversity(pick(starts_with("taxon")),index="shannon"),
        simpson = diversity(pick(starts_with("taxon")),index="simpson"),
        richness = specnumber(pick(starts_with("taxon"))),
        trueshannon = exp(shannon),
        truesimpson = 1/(1-simpson)
      ) %>%
      # get rid of reads columns
      select(-starts_with("Zotu"))
    
    # do plotting for each diversity metric
    metrics <- c(
      "Shannon" = "shannon",
      "True Shannon" = "trueshannon",
      "Simpson" = "simpson",
      "True Simpson" = "truesimpson",
      "zOTU Richness" = "richness"
    )
    metrics %>%
      imap(~{
        ggplot(div,aes(x=marker,y=.data[[.x]])) +
          geom_boxplot(aes(fill=marker)) +
          geom_point(aes(color = marker, fill = marker), shape = 21, size = 6, color = "black") +
          scale_fill_manual(values=fill_alpha(marker_pal,0.7)) + 
          xlab("Marker") + 
          ylab(.y) +
          theme_bw() +
          stat_summary(fun.data=data_summary) +
          theme(legend.position = "none",
                panel.grid = element_blank(),
                axis.title.y = element_text(size = 16),
                axis.title.x = element_text(size = 16),
                axis.text=element_text(size=14)) +
          coord_cartesian(ylim=c(0,max(div %>% pull(.x)))) + 
          labs(title=dsn)
      }) %>%
      # clean up the names of the returned list
      clean_names()
  })

# plot true shannon across the three datasets
div_plotz$sharkpen$true_shannon + div_plotz$aquarium$true_shannon + div_plotz$mock$true_shannon + plot_layout(axis_titles = "collect")

# zotu intersections (upset/venn) -----------------------------------------
source("upset.R")

# intersections
# which column to show intersections for
# for reasons, these have to be expressions, rather than strings
cols <- c(expr(zotu),expr(family))
# label maps
ll <- c("zotu" = "zOTUs", "family" = "families") 

upset_plotz <- 
  datasets %>%
  map(~{
    ds <- .x
    
    cols %>%
      set_names() %>%
      imap(~{
        ds %>%
          # arrange(desc(marker)) %>%
          pivot_wider(id_cols={{.x}},names_from="marker",values_from="reads",values_fill = 0,values_fn=~as.integer(sum(.x) > 0)) %>%
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
          labs(title=.x)
      })
  })

# show upset plots for family and zotu in the mock community
upset_plotz$mock$zotu / upset_plotz$mock$family


# expected vs unexpected species detections -------------------------------

# get ncbi ranked lineage dump
lineage <- get_lineage()

# create
expected_plotz <- datasets %>%
  imap(~{
    ds <- .x
    dsn <- .y
    if (file_exists(path(data_dir,str_glue("{dsn}_species.csv")))) {
      # read expected species from file and join in lineage info from ncbi taxonomy
      spp <- read_csv(path(data_dir,str_glue("{dsn}_species.csv")),col_types=cols()) %>%
        left_join(lineage,by=c("species" = "taxon")) %>%
        select(domain:genus,species) %>%
        arrange(across(domain:species))
      
      # get rid of unidentified things and collapse by species
      sp_id <- ds %>%
        filter(species != "unidentified") %>%
        group_by(marker,across(domain:species)) %>%
        summarise(reads = sum(reads)) %>%
        ungroup()   
      
      # create list of marker names to use with bind_cols below
      marker_cols <- markers %>%
        set_names() %>%
        map(~NA)
      
      # list of taxonomy names for bind_cols below
        tax_cols <- spp %>%
          select(domain:genus) %>%
          names() %>%
          set_names() %>%
          map(~NA_character_)
        
      # join everything together into the expected/unexpected table
      ss <- spp %>%
        # first set all known species to expected
        mutate(expected = TRUE) %>%
        # join in all the species we actually detected
        full_join(sp_id %>% distinct(across(domain:species)),by="species") %>%
        # the join will produce two versions of each level  down to genus (e.g., domain.x and domain.y)
        # here we add back the columns without suffixes
        bind_cols(tax_cols) %>%
        # now we set those taxonomy columns to whichever version (.x or .y) isn't NA
        # `coalesce` is the function that does this, and `coacross` is a wrapper that allows
        # us to use it in `across`
        mutate( across( all_of(names(tax_cols)), ~coacross(starts_with(cur_column())) ) ) %>%
        # get rid of the .x and .y variants of each taxonomic level
        select(-ends_with(".x"),-ends_with(".y")) %>%
        # any NA expected value is an unexpected species
        mutate(expected = replace_na(expected,FALSE)) %>%
        # add individual marker columns
        bind_cols(marker_cols) %>%
        # now get detections for each individual marker
        mutate( across( all_of(markers), ~species %in% (sp_id %>% filter(marker == cur_column()) %>% pull(species)) ) ) %>%
        # select columns into correct order
        select(domain,kingdom,phylum,class,order,family,genus,species,expected,everything()) %>%
        # sort by expected and taxonomy
        arrange(desc(expected),across(family:species)) %>%
        # add row number to sort species names
        mutate(n=row_number()) %>%
        # order species names as factor levels so they display in the order we want
        mutate(species = fct_reorder(species,-n)) %>%
        # get rid of temp column
        select(-n)  %>%
        # switch to long format
        pivot_longer(-c(domain:expected),names_to="marker",values_to="present")  %>%
        # create columns for plotting
        mutate( size = if_else(expected,"Expected species","Unexpected species") ) %>%
        # smash together marker and present true/false
        unite("clr",marker,present,remove = FALSE)
      
      # make a color palette so that each present marker gets its own
      # color, but all absent markers are just white
      mp <- c(
        marker_pal %>% set_names(str_c(names(marker_pal),"_TRUE")),
        rep("white",6) %>% set_names(str_c(markers,"_FALSE"))
      )
      
      ggplot(ss) + 
        geom_point(aes(x=marker,y=species,fill=clr,size=size),shape=21,color="black") +
        scale_size_discrete(guide="none") + 
        scale_fill_manual(values=mp,guide="none",na.value = "white") + 
        scale_x_discrete(position="top") +
        facet_wrap(~size,scales = "free") +
        theme_bw() + 
        theme(
          strip.background = element_rect(fill="#eeeeee"),
          strip.text = element_text(face = "bold")
        ) + 
        labs(title=dsn)
    }
  })

expected_plotz$mock
expected_plotz$aquarium


# extra garbage -----------------------------------------------------------
