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
save_pdf <- TRUE

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
  mock = "_mc[0-9]+$",
  aquarium = "_wa[0-9]+$",
  lagoon = "_sp[0-9]+$"
)

# map dataset names to display titles
title_map <- c(
  mock = "Mock community",
  mc = "Mock community",
  aquarium = "Waikīkī Aquarium",
  wa = "Waikīkī Aquarium",
  lagoon = "Coral Reef Lagoon",
  sp = "Coral Reef Lagoon"
)

# marker name map
marker_map <- c(
  "Berry 16S",
  "Leray CO1",
  "MiFish_E 12S",
  "MiFish_U 12S",
  "Riaz 12S",
  "Any marker"
) %>%
  set_names(c(markers,"Any marker"))

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

# helper to make nice number formats
num <- function(x,a=NULL) scales::label_comma(accuracy=a)(x)

# helper to relevel taxonomic factors
relevel_unid <- function(f) {
  unid <- grep("unidentified",f,value=TRUE) %>%
    unique() %>%
    sort()
  fct_relevel(f,unid,after=Inf)
}

# shortcut to do sum(x,na.rm=TRUE)
nasum <- function(x) sum(x,na.rm=TRUE)

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
      summarise(
        # sum sample reads
        across(all_of(samples),nasum),
        # pull out representative zotu (first one), and specify which marker it is
        representative = str_glue("{marker}:{otu[1]}"),
        
        # concatenate zotus with marker name like this marker(zotu1,zotu2,zotu3,...)
        zotus = str_glue("{marker}({str_c(otu,collapse=',')})"), 
        zotu_list = list(otu),
        zotu_count = n(),
        marker = marker
      ) %>%
      ungroup() %>%
      select(domain:species,marker,representative,zotus,zotu_list,zotu_count,all_of(samples))
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
    zotu_list = list(zotu_list),
    # sum read counts ()
    across(all_of(all_samples),nasum)
  ) %>%
  ungroup() %>%
  mutate(zotu_list = map(zotu_list,~set_names(.x,markers)))

# pull out sample names
samples <- all_samples

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
rn <- {if (rarefy) c("complete","rarefied") else c("complete")}

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

all_taxa <- datasets$complete %>%
  map(~.x %>% select(domain:species)) %>%
  list_rbind() %>%
  distinct(domain,kingdom,phylum,class,order,family,genus,species,.keep_all = TRUE) %>%
  arrange(class,family,species) 

# make a consistent color palette for markers
marker_pal <- paletteer_d("ggthemes::Tableau_10",n=length(markers)) %>%
  as.character() %>%
  append("black") 
marker_pal <- marker_pal %>%
  set_names(c(markers,"Any marker"))

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

# start here --------------------------------------------------------------

# this is just a placeholder space that allows us to jump here using rstudio's navigation shortcut

# generate tables ---------------------------------------------------------

# get raw sequencing data, including negative control filtering
raw_seq_data <- markers %>%
  set_names() %>%
  map(~{
    # save marker name into something nicer than `.x`
    marker <- .x
    
    # read marker data and filter to just chordates
    read_tsv(path(data_dir,str_glue("{marker}_data.tsv")),col_types = cols()) %>%
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
      select(-starts_with(marker),contains("blank"),where(~is.numeric(.x) && sum(.x) > 0)) %>%
      pivot_longer(-c(domain:seq_length),names_to = "sample",values_to="reads") %>%
      # make categorical groupings for bar plots
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
      # give them a sample type
      mutate(
        sample_type = case_when(
          sample == "blanks" ~ "Blanks",
          str_detect(sample,regex(dataset_map$lagoon,ignore_case = TRUE)) ~ "Coral Reef Lagoon",
          str_detect(sample,regex(dataset_map$aquarium,ignore_case = TRUE)) ~ "Waikīkī Aquarium",
          str_detect(sample,regex(dataset_map$mock,ignore_case = TRUE)) ~ "Mock community"
        )
      ) %>%
      # now get reads and zotus
      group_by(sample_type,sample,grouping) %>%
      summarise(reads = sum(reads), zotus = n_distinct(OTU)) %>%
      ungroup() %>%
      # get relative reads and relative zotus
      group_by(sample) %>%
      mutate(rel_reads = reads / sum(reads), rel_zotus = zotus/sum(zotus)) %>%
      # relevel sample types into order they are discussed
      mutate(sample_type = factor(sample_type, levels = c("Mock community", "Waikīkī Aquarium", "Coral Reef Lagoon"))) %>%
      arrange(sample_type,grouping)
  })

# make summary tables for each marker
marker_summary_tables <- raw_seq_data %>%
  map(~{
    .x %>%
      # get rid of blanks
      filter(sample_type != "Blanks") %>%
      # combine sharks & bony fishes into just "fishes"
      mutate(grp = fct_other(grouping,drop=c("Sharks & rays","Bony fishes"),other_level = "fishes")) %>%
      # get totals for each sample for all/fish reads/zotus
      group_by(sample_type,sample) %>%
      summarise(
        fish = sum(reads[grp == "fishes"]),
        reads = sum(reads),
        fish_zotus = sum(zotus[grp == "fishes"]),
        zotus = sum(zotus),
      ) %>%
      ungroup() %>%
      # now further collapse by sample type
      # and get nice-looking summary stats
      group_by(sample_type) %>%
      summarise(
        `Total reads` = num(sum(reads)),
        `Mean reads per sample` = str_glue("{num(mean(reads))} ± {num(sd(reads))}"),
        
        `Total fish reads` = num(sum(fish)),
        `Mean fish reads per sample` = str_glue("{num(mean(fish))} ± {num(sd(fish))}"),
        
        `Mean ZOTUs per sample` = str_glue("{num(mean(zotus))} ± {num(sd(zotus))}"),
        
        `Mean fish ZOTUs per sample` = str_glue("{num(mean(fish_zotus))} ± {num(sd(fish_zotus))}"),
      ) %>%
      ungroup() %>%
      rename(`Sample type` = sample_type)
  })

# save summary tables
marker_summary_tables %>%
  iwalk(~{
    write_tsv(.x,path(tbl_dir,str_glue("{.y}_marker_summary.tsv")))  
  })

# hacky function to make things like "family" into "families" and "genus" into "genera"
pluralize <- function(s) {
  map_chr(s,\(s) {
    if (str_detect(s,regex("y$",ignore_case = TRUE))) {
      str_replace(s,"y$","ies")
    } else if (str_detect(s,regex("um$",ignore_case = TRUE))) {
      str_replace(s,"um$","a")
    } else if (s != "ZOTUs" & str_detect(s,regex("us$",ignore_case = TRUE))) {
      str_replace(s,"us$","era")
    } else if (s == "class") {
      return("classes")
    } else if (s %in% c("species","ZOTUs")) {
      return(s)
    } else {
      str_c(s,"s")
    }
  })
}


reverse_marker_map <- markers %>%
  set_names(make_clean_names(markers))

all_summary  <- fishes_filtered %>%
  # pivot sample names to long with marker, type, and sample #
  pivot_longer(any_of(samples),names_pattern = '^(.+)_(..)([^_]+)$',names_to = c("marker","type","sample"),values_to = "reads") %>%
  # keep only detected entries
  filter(reads > 0) %>%
  mutate(
    # rename marker to something nice
    marker = reverse_marker_map[marker],
    # filter to the zotus of the marker
    zotu_list = map2(zotu_list,marker,\(zl,m) zl[[m]]),
    # rename sample type to something nice
    type = title_map[type],
    # de-factorize family to species
    across(family:species,as.character)
  ) %>%
  # unlist the zotu list
  unnest(zotu_list) %>%
  # get per-sample numbers for each marker and sample type
  group_by(marker,type,sample) %>%
  summarise(
    # zotus per replicate
    ZOTUs_per_sample = n_distinct(zotu_list),
    # re-list zotus
    zotu_list = list(zotu_list),
    # families, etc. per replicate
    across(family:species,~n_distinct(identified(.x)),.names = "{.col}_per_sample"),
    # re-list families, etc.
    across(family:species,list)
  ) %>%
  # get number per marker and sample type
  group_by(marker,type) %>%
  summarise(
    # get mean per-sample families, zotus, etc
    across(ends_with("per_sample"),mean,.names = "{.col}_mean"),
    # get per-sample families, zotus, etc std. dev.
    across(ends_with("per_sample"),sd,.names = "{.col}_sd"),
    # overall unique zotus
    ZOTUs = n_distinct(flatten(zotu_list)),
    # overall unique identified families, etc.
    across(family:species,~n_distinct(identified(flatten(.x))))
  ) %>%
  # sort by sample type and marker
  rename_with(.cols=ZOTUs:species,pluralize) %>%
  pivot_longer(contains("per_sample"),names_pattern = "(.+)_per_sample_(.+)", names_to = c("taxon","measurement")) %>%
  mutate(value = round(value,1)) %>%
  group_by(marker,taxon,type) %>%
  summarise(
    across(-c(value,measurement),first),
    value = str_c(value,collapse = " ± ")
  ) %>%
  mutate(taxon = str_c(pluralize(taxon)," per replicate")) %>%
  pivot_wider(names_from="taxon",values_from="value") %>%
  rename_with(.cols = matches("^(fam|gen|spe)"),str_to_sentence) %>%
  mutate(type = factor(type,levels=unique(title_map))) %>%
  arrange(type,marker) %>%
  select(`Sample type`=type,Marker=marker,starts_with("fam"),starts_with("gen"),starts_with("spe"),starts_with("zotu"))


write_tsv(all_summary,path(tbl_dir,"all_summary.tsv"))

# community multivariate statistics and pcoa plots (Fig. pt. 1) ------------

community_stats <- datasets %>%
  map(~{
    dataset <- .x
    dataset %>%
      imap(~{
        dsn <- .y
        div <- .x %>%
          { if (filter_unid) filter_unidentified(., "species") else . } %>%
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
          decostand("pa")
        
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

## Save permanova summary results
community_stats %>%
  iwalk(~{
    r <- .y
    .x %>%
      iwalk(~{
        p <- .y
        .x$permanova %>% rownames_to_column(var = "x") %>%
        write_tsv(path(tbl_dir,str_glue("permanova_{p}_{r}.tsv")))
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
        sd <- datasets %>%
          pluck(ds) %>%
          pluck(dsn) %>%
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
          left_join(sd,by="sample") %>%
          mutate(marker = factor(marker_map[marker],levels=marker_map))
        
        # plot pcoa
        ggplot(bb,aes(x=PCoA1,y=PCoA2,fill=marker)) + 
          geom_point(shape=21,size=4,color="black",show.legend=TRUE) + 
          scale_fill_manual(values=marker_pal %>% set_names(marker_map),name="Primer",drop=FALSE) +
          theme_bw() +
          xlab(str_glue("PCoA1 ({scales::percent(axes[1],accuracy=0.1)})")) + 
          ylab(str_glue("PCoA2 ({scales::percent(axes[2],accuracy=0.1)})"))
      })
  })

# chord plots (Fig. 1 pt. 2) ----------------------------------------------
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

chord_plotz <- datasets$complete %>%
  imap(~{
    dsn <- .y
    dd <- .x %>%
      { if (filter_unid) filter_unidentified(., "family") else . } %>%
      select(-sample, -type) %>% distinct() %>%
      arrange(class, family, marker) %>%
      select(marker, family, marker_zotu_count)
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
    as.ggplot(function() make_chord(dd,c(1e-9,1,1e-9,1e-9),c(marker_pal,palettes$family),units="in",group=grp))
  })

comp_layout = c(
  patchwork::area(l = 72, r = 96, t = 5, b = 7),
  patchwork::area(l = 72, r = 96, t = 13, b = 15),
  patchwork::area(l = 72, r = 96, t = 21, b = 23),
  patchwork::area(l = 1, r = 72, t = 1, b = 10),
  patchwork::area(l = 1, r = 72, t = 9, b = 18),
  patchwork::area(l = 1, r = 72, t = 17, b = 26)
  )

tag_positions = list(c(0,1), c(0,1), c(0,1),
                     c(0.075,0.8), c(0.075,0.8), c(0.075,0.8))

chord_pcoa_composites <- pcoa_plotz %>%
  map(~{
    plots = map2(c(.x, chord_plotz), tag_positions, ~ .x + theme(legend.position = "none",
                                                                 plot.tag.position = .y))
    wrap_plots(plots, ncol = 2, byrow = F) +
      plot_layout(design = comp_layout) +
      plot_annotation(theme = theme(plot.margin = margin(t = -75, r = 20, b = -25)),
                      tag_levels = list(c("(d)", "(e)", "(f)",
                                          "(a)", "(b)", "(c)")))
  })

if (save_pdf) {
  chord_pcoa_composites %>%
    iwalk(~ggsave(path(fig_dir,str_glue("chord_pcoa_composite_{.y}.pdf")),.x,device=cairo_pdf,width=12,height=22,units="in"))
}

# expected vs unexpected species detections (Figs. 2-4) -------------------

# load known invasive/introduced species
ais <- read_csv(path(data_dir,"ais.csv"))

# hacky little config map to say whether we want to see ALL possible expected species in the figure (or not)
# this is done so that the lagoon figure doesn't show ~1200 speces in the 'expected' column
show_all <- c(
  mock = TRUE,
  aquarium = TRUE,
  lagoon = FALSE
)

# hacky little config map to say whether we want to see invasive/introduced species
# (of course we don't care about this for the aquarium dataset)
show_inv <- c(
  mock = FALSE,
  aquarium = FALSE,
  lagoon = TRUE
)

exp_map = list(
  lagoon = c("Not known from Hawai‘i","Known from Hawai‘i"),
  aquarium = c("Not known from aquarium tank","In aquarium tank"),
  mock = c("Not in mock community","In mock community")
)

# whether to color invasive/introduced/whatever
color_introduced <- TRUE

# whether to do point or tile
tile <- TRUE

ranks <- c("genus","species")

# create
expected_plotz <- datasets %>%
  map(~{
    ranks %>%
      set_names() %>%
      map(\(taxon) {
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
                filter_unidentified(taxon) %>%
                # group by marker and species
                group_by(marker,across(domain:all_of(taxon))) %>%
                summarise(reads = sum(reads)) %>%
                ungroup() %>%
                # keep only family to species
                select(-c(domain:order),class) %>%
                # switch to wide (marker x detected)
                pivot_wider(names_from="marker",values_from=reads,values_fn=~.x > 0,values_fill=FALSE) 
              
              # glue text for invasive species/genera
              inv_spp <- str_glue('{{{taxon}}}<sup>**†**</sup>')
              nat_gen <- str_glue('{{{taxon}}}<sup>**\\***</sup>')
              
              # get join columns
              join_cols <- sp_id %>%
                select(-any_of(markers)) %>%
                names()
              
              tax_cols <- ds %>%
                select(-c(marker:last_col())) %>%
                names()
              
              join_cols <- tax_cols[sort(match(join_cols,tax_cols))]
              
              # join everything together into the expected/unexpected table
              ss <- spp %>%
                # keep only columns we want
                select(all_of(join_cols)) %>%
                distinct(across(all_of(join_cols))) %>%
                # first set all expected species to expected
                mutate(expected = TRUE) %>% 
                # join in species we actually detected
                full_join(sp_id,by=join_cols) %>%
                select(all_of(join_cols),expected,any_of(markers)) %>%
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
                pivot_longer(-c(all_of(join_cols),expected),names_to="marker",values_to="present")  %>%
                # drop markers with zero detections
                group_by(marker) %>%
                filter(sum(present) > 0) %>%
                ungroup() %>%
                # filter out unexpected 'Any' markers
                filter(!(!expected & marker == "Any marker")) %>%
                # create columns for plotting
                mutate( 
                  exp = exp_map[[dsn]][as.integer(expected)+1]
                ) %>%
                # smash together marker and present true/false
                mutate(marker = marker_map[marker]) %>%
                unite("clr",marker,present,remove = FALSE) %>%
                { if(show_all[dsn]) mutate(.,marker = fct_relevel(marker,"Any marker",after=Inf)) else . } %>%
                # annotate known introduced/invasive species
                mutate(
                  "{taxon}" := case_when(
                    color_introduced &
                      taxon == "species" & 
                      show_inv[dsn] &
                      .data[[taxon]] %in% ais[[taxon]] ~
                      str_glue(inv_spp),
                    color_introduced &
                      taxon == "species" & 
                      show_inv[dsn] &
                      !expected &
                      !(.data[[taxon]] %in% ais[[taxon]]) &
                      genus %in% ais$genus &
                      genus %in% spp$genus ~ str_glue(nat_gen),
                    .default = .data[[taxon]]
                  ),
                  "{taxon}" := str_glue(str_glue("*{{{taxon}}}*"))
                ) %>%
                # sort by expected and taxonomy
                arrange(desc(expected),class,across(family:all_of(taxon))) %>%
                # add row number to sort species names
                mutate(n=row_number()) %>%
                # order species names as factor levels so they display in the order we want
                mutate("{taxon}" := fct_reorder(.data[[taxon]],-n)) %>%
                # get rid of temp column
                select(-n)  
              
              # make a color palette so that each present marker gets its own
              # color, but all absent markers are just white
              mp = c(
                marker_pal %>% set_names(str_c(marker_map, "_TRUE")),
                rep("white",length(marker_pal)) %>% set_names(str_c(marker_map, "_FALSE"))
              )
              
              # do plot
              ggplot(ss) + 
                # points/tiles for presence/absence of species x marker
                { if (tile) geom_tile(aes(x=marker,y=.data[[taxon]],fill=clr),color="gray10") else geom_point(aes(x=marker_map[marker],y=.data[[taxon]],fill=clr),shape=21,color="black",size=2.5) } +
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
      }) # end ranks
  })

# show a couple example plots
# expected_plotz$complete$mock
# expected_plotz$complete$aquarium
# expected_plotz$complete$lagoon

# another hacky little map so we get the right plot sizes
plotsizes <- list(
  mock = c(9.8,13),
  aquarium = c(7.7,7.7),
  lagoon = c(7.7,7.7)
)

# save expected plot figures
if (save_pdf) {
  expected_plotz %>%
    iwalk(~{
      r <- .y
      .x %>% iwalk(~{
        t <- .y
        .x %>% iwalk(~{
          p <- .y
          ggsave(path(fig_dir,str_glue("expected_{p}_{t}_{r}.pdf")),.x,device=cairo_pdf,width=plotsizes[[p]][1],height=plotsizes[[p]][2],units="in")
        })
      })
    })
}

# primer performance plots (Fig. 5) ---------------------------------------
primer_taxa_pal <- c(
  "Bony fishes" = "#3B01FE",
  "Sharks & rays" = "#A9A9A9",
  "Mammals" = "#EE82EF",
  "Bacteria" = "#ED0105",
  "Other" = "#F2A503"
)

primer_plotz <- raw_seq_data %>%
  imap(~{
    pd <- .x %>%
      # get rid of blanks
      filter(sample_type != "Blanks")
    nn = .y
    
    read_rel <- ggplot(pd) + 
      geom_col(aes(x=sample,y=rel_reads,fill=grouping), show.legend = TRUE) + 
      scale_fill_manual(values = primer_taxa_pal, name = "Group", drop=FALSE) + 
      scale_y_continuous(labels=scales::percent, expand = expansion(mult=c(0.01,NA))) +
      labs(subtitle = marker_map[.y], y="Relative sequence abundance") + 
      facet_wrap(~sample_type,scales="free_x",strip.position = "bottom") &
      theme_bw() & 
      theme(
        panel.border = element_blank(),
        axis.line = element_line(color="black"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(face="bold"),
        strip.background = element_rect(fill=NA,color=NA)
      )
  })

primers_rel_reads_plot <- primer_plotz %>%
  map(~.x) %>%
  reduce(`+`) + 
  plot_layout(ncol = 2, axis_titles = "collect", guides = "collect") +
  plot_annotation(tag_levels = plot_tags) & 
  labs(x="")

if (save_pdf) {
  ggsave(path(fig_dir,"primer_relative_reads.pdf"),device=cairo_pdf,primers_rel_reads_plot,width=9.5,height=6.5,units="in")
}

# family intersections (upset plots) (Fig. 6) -------------------------------

upset_plotz <- datasets %>%
  map(~{
    .x %>%
      map(~{
        ds <- .x
        ds %>%
          filter_unidentified(., "family") %>%
          # create a presence/absence dataset where we pivot to wider by marker x family
          # and summarise each occurrence as a zero or one (1 = sum(reads) > 0)
          pivot_wider(id_cols=family, names_from="marker", values_from="reads", values_fill = 0, values_fn=~as.integer(sum(.x) > 0)) %>%
          # plot the upset plot
          upset_plot(
            name_column=family,
            data_columns=any_of(markers),
            label_top_bars = TRUE,
            label_side_bars = TRUE,
            bar_lab = str_glue("Families"),
            sidebar_lab = "Families", 
            group_palette = marker_pal
          ) %>%
          wrap_elements()
      })
  })

# show upset plots for family and zotu in the unrarefied mock community
# upset_plotz$complete$lagoon$family / upset_plotz$complete$lagoon$zotu + plot_annotation(tag_levels = plot_tags)
# upset_plotz$complete$aquarium$family / upset_plotz$complete$aquarium$zotu + plot_annotation(tag_levels = plot_tags)
# upset_plotz$complete$mock$family / upset_plotz$complete$mock$zotu + plot_annotation(tag_levels = plot_tags)

# put the plots together for each sample type
upset_composites <- upset_plotz %>%
  map(~{
    reduce(.x, `/`) + plot_annotation(tag_levels = plot_tags)
  })

# save them
if (save_pdf) {
  upset_composites %>%
    iwalk(~{
      r <- .y
      ggsave(
        filename = file.path(fig_dir, str_glue("upset_{r}.pdf")),
        plot = .x,
        device = cairo_pdf,
        width = 10,
        height = 13,
        units = "in"
      )
    })
}

# taxon bar plots (Figs. S1-S2) -------------------------------------------

# taxonomic levels at which to group bar plots
plot_levels <- c("order","family")

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
              { if (pl %in% c("class", "order")) . else filter_unidentified(., pl) } %>%
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
              geom_col(aes(x=marker_map[marker],y=rel,fill=.data[[pl]]),show.legend = TRUE) +
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

# rel_taxon_plotz$complete$family$aquarium
# rel_taxon_plotz$complete$family$mock
# rel_taxon_plotz$complete$family$lagoon

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
# rel_taxon_composites$complete$family
# rel_taxon_composites$complete$order

# relative zotu abundance heatmaps (Figs. S3-S4) ---------------------------

# plot zotus by taxonomic level
plot_levels <- c("order","family")

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
              { if (pl %in% c("class", "order")) . else filter_unidentified(., pl) } %>%
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
              mutate("{pl}" := fct_relevel(.data[[pl]],"unidentified")) %>%
              suppressWarnings()
            
            ggplot(dd) + 
              geom_tile(aes(x=marker_map[marker],y=.data[[pl]],fill=rel),color="grey8") + 
              scale_fill_paletteer_c("viridis::turbo",name="Relative\nabundance")  +
              scale_y_discrete(limits=rev) + 
              labs(x="Marker",y=str_to_sentence(pl)) + 
              theme(axis.title = element_text(face="bold"))
          })
      })
  })

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


# proportion of expected taxa (not included in MS) ------------------------

proportion_expected_plotz <- datasets %>%
  map(~{
    plot_data <- .x %>%
      # dump the lagoon since we have something like 1,200 "expected" species
      discard_at("lagoon") %>%
      imap(~{
        ds <- .x
        dsn <- .y
        if (file_exists(path(data_dir,str_glue("{dsn}_species.csv")))) {
          # read expected species from file and join in lineage info from ncbi taxonomy
          # we use the lineage info to sort by family in addition to genus/species
          spp <- read_csv(path(data_dir,str_glue("{dsn}_species.csv")),col_types=cols()) %>%
            mutate(expected=TRUE) %>%
            select(family:species,expected)
          # prepare list of detected species
          sp_id <- ds %>%
            # drop unidentified things
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
          # join everything together into the expected/unexpected table
          plot_data <- spp %>%
            # first set all expected species to expected
            # join in species we actually detected
            full_join(sp_id,by=c("family","genus","species")) %>%
            mutate(
              # any NA expected value is an unexpected species
              expected = replace_na(expected,FALSE),
              # any NA marker value is a negative detection
              across(any_of(markers),~replace_na(.x,FALSE))
            )  %>%
            # switch to long format
            mutate(`Any marker` = rowSums(pick(any_of(markers))) > 0) %>%
            pivot_longer(c(any_of(markers),`Any marker`),names_to = "marker",values_to = "present") %>%
            group_by(marker) %>%
            summarise(
              present = sum(present & expected),
              expected = sum(expected)
            ) %>%
            ungroup() %>%
            mutate(
              proportion = present/expected,
              dataset=title_map[dsn],
              marker = marker_map[marker],
              marker = fct_relevel(marker,"Any marker",after=Inf)
            ) 
          
          ggplot(plot_data) +
            geom_col(aes(x=marker,y=proportion,fill=names(marker_map)[marker])) + 
            scale_y_continuous(limits=c(0,1),expand = c(0,NA), labels = scales::percent) + 
            scale_fill_manual(values=marker_pal,name="Primer",guide="none") + 
            theme_bw() +
            theme(
              panel.border = element_blank(),
              axis.line = element_line(color="black"),
              panel.grid = element_blank(),
              axis.text.x = element_text(angle=30, hjust=1),
              axis.title = element_text(face="bold")
            ) + 
            labs(x="Primer",y="% expected species recovered")
        }
      }) %>%
      reduce(`+`) + 
      plot_layout(guides = "collect", axis_titles = "collect") + 
      plot_annotation(tag_levels = plot_tags)
  })

if (save_pdf) {
  proportion_expected_plotz %>%
    iwalk(~ggsave(path(fig_dir,str_glue("proportion_expected_{.y}.pdf")),.x,device=cairo_pdf,width=9,height=4,units="in"))
}
