#' Plot the standardised and unstandardised indices
#' 
#' In this plot the unstandardised indices is the geometric mean of the data.
#' 
#' @param index A dataframe with unstan and stan idex for the model. 
#' @param year the year or time label.
#' @param fill the fill colour for the percentiles.
#' @param probs The percentiles to be computed by the \code{quantile} function.
#' @param rescale the index of the series to rescale to. If set to NULL then no rescaling is done.
#' @param show_unstandardised show the unstandardised series or not.
#' @return a \code{ggplot} object.
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @export
#' 
plot_index <- function(index, 
                       fill = "black", 
                       probs = c(0.25, 0.75),
                       rescale = 1,
                       predictor = NULL,
                       show_unstandardised = TRUE) {
  
  index_long <- index %>%
    pivot_longer(
      cols = -level,
      names_to = c("is_stan", "stat", "is_scaled"),
      names_pattern = "^(stan|unstan)(Lower|Upper)?_?(unscaled)?$"
    ) %>%
    mutate(
      stat = ifelse(stat == "" , "median", stat),
      is_scaled = ifelse(is_scaled == "unscaled", "Binomial - probabilities", "Relative index"),
      is_stan = factor(is_stan, levels = c("stan", "unstan"), 
                       labels = c("Standardised", "Unstandardised"))
      
    ) %>%
    # Pivot 'median', 'Lower', and 'Upper' back into their own columns
    pivot_wider(
      names_from = stat, 
      values_from = value
    )
  
  if (!show_unstandardised) {
    index_long <- filter(index_long, is_stan == "Standardised")
  }
  
  
  if (!("unstan_unscaled" %in% names(index))) {
    # No strip text if positive index
    facet_theme <- theme(strip.text = element_blank())
    # no unscaled component if positive index
    index_long <- filter(index_long, is_scaled == "Relative index")
  } else facet_theme <- theme(legend.position.inside = c(0.05, 0.35))
  
  
 p <- ggplot(index_long, aes(x = level, y = median, group = is_stan)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), 
                data = filter(index_long, is_stan=='Standardised'),
                fill = fill, color = NA, alpha = 0.1) +
    geom_line(aes(linetype = is_stan, colour = is_stan)) +
    geom_point(aes(shape = is_stan, colour = is_stan), size = 3) +
    facet_wrap(~is_scaled, ncol = 1, scales = "free_y")+
    labs(x = "Fishing year", y = "Index") +
    scale_shape_manual(values = c('Standardised' = 16, 'Unstandardised' = 1))+
    scale_colour_manual(values = c('Standardised' = fill, 'Unstandardised' = "grey40")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_cowplot() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "inside",
      legend.position.inside = c(0.05, 0.9),
      legend.direction = "horizontal",
      legend.justification = c(0, 0), 
      legend.title = element_blank(),
      legend.key.width = unit(2, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    ) +
    facet_theme 

return(p)
}


