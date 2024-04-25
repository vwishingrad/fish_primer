# I don't love the pairwise adonis options available in packages so I made my own
pairwise_adonis <- function(comm, factors, permutations = 1000, correction = "fdr", method = "bray") {
  # get possible pairwise factor combinations
  factor_combos <- combn(unique(as.character(factors)), 2)
  # map through factor combinations and run model for each pair,
  # dumping the output into a tibble
  model_output <- map_dfr(array_branch(factor_combos,2), ~{
    fact <- factors[factors %in% .x]
    
    # get our factor specific community or distance matrix
    if (inherits(comm,'dist')) {
      dd <- as.dist(as.matrix(comm)[factors %in% .x, factors %in% .x])
    } else {
      comm <- as(comm,'matrix')
      dd <- vegdist(comm[factors %in% .x,], method = method)
    }
    
    # run the permanova model and return the results in a list
    # to make the data frame rows
    model <- adonis2(dd ~ fact, permutations = permutations) %>%
      as_tibble()
    df <- model$Df[1]
    ss <- model$SumOfSqs[1]
    pseudo_f <- model$F[1]
    p_val <- model$`Pr(>F)`[1]
    
    list(
      "left" = .x[1],
      "right" = .x[2],
      "df" = df,
      "ss" = ss,
      "pseudo_f" = pseudo_f,
      "p_val" = p_val
    )
  })
  
  # return the results with adjusted p-values
  model_output %>%
    mutate(
      p_adj = p.adjust(p_val, method = correction)
    )
}


# make an upset plot
upset_plot <- function(dataset, name_column, data_columns, dot_size = 6, line_size=2,
                       prefilter=F, intersects=NA, min_intersects=0, bar_lab = "intersections",
                       sidebar_lab = "number in category", label_top_bars = FALSE, label_side_bars = FALSE,
                       group_palette = NULL) {
  if (prefilter) {
    dataset <- dataset %>%
      mutate(
        across({{data_columns}},~if_else(.x > 0,1,0,missing=0))
      )
  }
  
  # get category names
  sets <- dataset %>%
    select({{data_columns}}) %>%
    names()
  
  dataset <- dataset %>%
    tidyr::unite("code",{{data_columns}},sep="",remove = FALSE) 
  
  dataset_long <- dataset %>%
    # pivot_longer({{data_columns}},names_to = "metric_name", values_to = "metric")
    pivot_longer(all_of(sets),names_to = "metric_name", values_to = "metric")
  
  data1 <- dataset_long %>%
    group_by(code,metric_name,metric) %>%
    summarise(n = n_distinct({{name_column}})) %>%
    arrange(n) %>%
    ungroup()
  
  data2 <- dataset_long %>%
    group_by(metric_name) %>%
    summarise(n=sum(metric))
  
  x_breaks <- dataset_long %>%
    group_by(code) %>%
    summarise(n=n_distinct({{name_column}})) %>%
    arrange(desc(n)) %>%
    filter(n > min_intersects) %>%
    pull(code) %>%
    as.character()
  
  if (!is.na(intersects)) {
    x_breaks <- x_breaks[1:intersects]
  }
  
  intersects <- length(x_breaks)
  
  top_bars <-
    ggplot(data1, aes(x=fct_reorder(factor(code),-n), y=n)) +
    geom_col(fill="grey5", position="dodge") +
    scale_x_discrete(limits=x_breaks) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    theme_bw() +
    ylab(bar_lab) +
    theme(legend.position = "none", 
          # axis.title = element_blank(),
          axis.line.y.left = element_line(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(0,0,0,0,"cm"))
  
  if (!is.null(group_palette)) {
    fill_vals <- group_palette
    if (all(is.na(names(fill_vals)))) {
      names(fill_vals) <- sets
    }
  } else {
    fill_vals <- rep("grey5",length(sets))
    names(fill_vals) <- sets
  }
  
  side_bars <-
    ggplot(data2, aes(x=metric_name, y=n)) +
    geom_col(aes(fill=metric_name), position="dodge") +
    scale_fill_manual(values=fill_vals) + 
    scale_x_discrete(
      position="top",
      limits=rev(sets)
    ) +
    scale_y_reverse(
      labels=scales::format_format(big.mark = ",", decimal.mark = ".", scientific = FALSE, digits=0),
      expand = expansion(mult = c(0.6, 0))
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      axis.line.x = element_line(),
      panel.border = element_blank(),
      legend.position = "none",
      axis.title = element_blank(),
      axis.text.x = element_text(size=6),
      axis.text.y = element_text(size=12),
      panel.grid = element_blank(),
      axis.title.x = element_text(),
      plot.margin = margin(0,0,0,0,"cm")
    ) +
    ylab(sidebar_lab) 
  
  if (label_side_bars) { 
    side_bars <- side_bars + 
      geom_text(aes(label=n), position = position_dodge(0.9),  hjust=1.1, vjust=0.5) 
  }
  
  if (label_top_bars) {
    top_bars <- top_bars +  
      geom_text(aes(label=n), position = position_dodge(0.9), hjust=0.5, vjust=-0.25)
  }

  dot_lines <- data1 %>%
    filter(metric == 1) %>%
    group_by(code) %>%
    summarise(
      f = list(factor(metric_name,levels=sets)),
      n = unique(n)
    ) %>%
    mutate(
      start = map_chr(f,~sets[min(as.integer(.x))]),
      end = map_chr(f,~sets[max(as.integer(.x))])
    ) %>%
    select(-f)
  
  cols <- c("0" = "grey77", "1" = "grey2")  
  data1 <- data1 %>%
      mutate(color_group = str_glue("{metric_name}_{metric}"))
  
  cols <- data1 %>%
    ungroup() %>%
    mutate(
      color_value = fill_vals[metric_name],
      color_value = case_when(
        metric == 0 ~ "grey95",
        TRUE ~ color_value
      )
    ) %>%
    distinct(color_group,color_value) %>%
    deframe()
  
  dots <-
    ggplot(data1, aes(y=metric_name, x=fct_reorder(factor(code),-n))) +
    geom_point(shape=21, size=dot_size, colour="black", aes(fill=color_group)) +
    geom_point(data=data1 %>% filter(metric == 1),shape=19,size=dot_size/2,color="black") + 
    geom_segment(data=dot_lines,mapping=aes(x=fct_reorder(factor(code),-n),xend=fct_reorder(factor(code),-n),y=start,yend=end),linewidth=line_size) +
    scale_fill_manual(values = cols) +
    scale_x_discrete(limits=x_breaks) +
    scale_y_discrete(limits=rev(sets)) +
    theme_minimal() +
    labs(x="", y="") +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.margin = margin(0,0,0,0,"cm"))
  layout <- "
    #1111
    #1111
    #1111
    23333
    23333
  "
  return(top_bars + side_bars + dots + plot_layout(design=layout))
}
