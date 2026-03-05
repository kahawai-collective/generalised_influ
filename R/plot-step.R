#' A Bayesian version of the step-plot
#' 
#' This requires that all steps be run using brms and then provided as a list of model fits.
#' 
#' @param fits a list of model fits in the order that you want to compare them
#' @param year the year or time label
#' @param fill the colour of the credible interval ribbon
#' @param probs the quantiles to plot
#' @param show_probs plot the quantiles or not
#' @import brms
#' @import ggplot2
#' @import dplyr
#' @export
#' 
plot_step <- function(fits, year = NULL, fill = "purple",
                      probs = c(0.25, 0.75), show_probs = TRUE) {
  
  m <- length(fits)
  
  fout <- list()
  for (i in 1:m) {
    fout[[i]] <- get_index(fit = fits[[i]], year = year, probs = probs)
  }
  
  df <- NULL
  df_dash <- NULL
  df_grey <- NULL
  for (i in 1:m) {
    df <- rbind(df, fout[[i]])
    if (i > 1) {
      xx <- fout[[i - 1]] %>% mutate(Model = fout[[i]]$Model, line = i)
      df_dash <- rbind(df_dash, xx)
    }
    if (i > 2) {
      xx <- fout[[i - 2]] %>% mutate(Model = fout[[i]]$Model, line = i)
      df_grey <- rbind(df_grey, xx) # bug - these needs to include all prevous models in grey (i.e. if 4 models are provided)
    }
  }
  
  p <- ggplot(data = df) +
    geom_line(data = df_grey, aes(x = .data$Year, y = .data$Median, group = .data$Model), colour = "grey", linetype = "solid") +
    geom_line(data = df_dash, aes(x = .data$Year, y = .data$Median, group = 1), colour = "black", linetype = "dashed")
  
  if (show_probs) {
    p <- p + geom_ribbon(data = df, aes(x = .data$Year, ymin = .data$Qlower, ymax = .data$Qupper, group = 1), alpha = 0.3, colour = NA, fill = fill)
  }
  
  p <- p + 
    geom_line(data = df, aes(x = .data$Year, y = .data$Median, group = 1)) +
    geom_point(data = df, aes(x = .data$Year, y = .data$Median)) +
    facet_wrap(Model ~ ., ncol = 1, strip.position = "top") +
    labs(x = NULL, y = "Index") +
    scale_fill_manual(values = c(NA, "black")) +
    scale_colour_manual(values = c("black", "black")) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_shape_manual(values = c(NA, 19)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing.y = unit(0, "lines"))
  
  return(p)
}

#' This requires that all steps be run using brms and then provided as a list of model fits.
#' 
#' @param step_df dataframe generated weith get_step() function
#' @import ggplot2
#' @import dplyr
#' @export
#' 
#' 
#' 

plot_step2 <- function(step_df, compare_step_df = NULL){
  
  # Identify position of first model
  
  start_idx <- match("stanUpper", names(step_df)) +1
  
  model_names <- names(step_df)[start_idx:ncol(step_df)]
  
  df_long <- step_df %>%
    dplyr::select(Year = level, all_of(model_names)) %>%
    pivot_longer(-Year, names_to = "Model", values_to = "Index") %>%
    mutate(Model = factor(Model, levels = model_names)) # Keeps them in order
  
  # Generate the background data for every facet
  # For each step, we take all data from that step and all previous steps
  df_all_steps <- lapply(seq_along(model_names), function(i) {
    df_long %>%
      filter(as.numeric(Model) <= i) %>%
      mutate(
        FacetTarget = factor(model_names[i]),
        LineType = case_when(
          as.numeric(Model) == i     ~ "Current",
          as.numeric(Model) == i - 1 ~ "Previous",
          TRUE                       ~ "Historical"
        )
      )
  })%>%bind_rows()
  
  # Comparison logic
  if(!is.null(compare_step_df)){
    
    df_all_steps<- df_all_steps %>%
      bind_rows(compare_step_df %>%
                  dplyr::select(Year = level, all_of(model_names)) %>%
                  pivot_longer(-Year, names_to = "Model", values_to = "Index") %>%
                  mutate(Model = factor(Model, levels = model_names),                        # Keeps them in order
                         FacetTarget = Model,
                         LineType = 'Compare'))
    
  }
  # End of comparison logic
  
  p <- ggplot(df_all_steps, aes(x = Year, y = Index, group = interaction(Model, LineType))) +
    geom_line(aes(color = LineType, linetype = LineType == "Previous"), linewidth = 0.5 ) +
    scale_color_manual(values = c("Historical" = "grey85", "Previous" = "black", "Current" = "royalblue", "Compare" = 'red')) +
    geom_point(data = filter(df_all_steps, LineType == "Current"), 
               color = "royalblue") +
    
    
    # geom_label(data = df_all_steps %>%
    #              group_by(FacetTarget) %>%
    #              summarize(Model = last(Model)), 
    #            aes(label = FacetTarget, x = -Inf, y = -Inf),
    #            # Fixed alignment settings:
    #            hjust = 0,             
    #            vjust = 0,             
    #            label.padding = unit(1, "lines"), 
    #            label.size = NA,       
    #            fill = NA,            
    #            size = 3.5) +
    labs(x = "Fishing year", y = "Index") +
    scale_y_continuous(limits = c(0.6, NA), 
                       expand = expansion(mult = c(0, 0.1)), 
                       guide = guide_axis(check.overlap = TRUE),
                       breaks = function(x) unique(pretty(x)[pretty(x) != 0]) )+
    facet_wrap(~FacetTarget, ncol = 1) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      strip.background = element_blank(), 
      strip.text = element_blank(),       
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
      panel.spacing = unit(0, "lines"),
      axis.text = element_text(size = 10)
    )
  
  return(p)
}
