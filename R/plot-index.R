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
                       component = 'Combined',
                       fill = "black", 
                       probs = c(0.25, 0.75),
                       rescale = 1,
                       show_unstandardised = TRUE,
                       overlay = FALSE) {
  
  # Default component "Combined", overwrite it if it does not exist
  if(length(unique(index$Index))==1) component <- unique(index$Index)
  
  # Filter to only one index (bin or pos) if we are not looking for combined
  if (component!= "Combined") {
    index <- index %>%
      filter(Index == component)
  }
  
  # Conditionally remove unstandardised idx
  if (!show_unstandardised|component=='Combined') {
    index <- filter(index, is_stan)
  }
  
  #no unscaled component if not binomial index
  if (component!='Binomial') {
    index <- filter(index, is_scaled)
  }
  
  if (overlay){
    cols <- c('dodgerblue', '#F5B915FF', '#08235FFF' )
    p <- ggplot(index, aes(x = level, y = median, group = Index, colour = Index, fill = Index)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper),
                  data = filter(index, is_stan),
                  color = NA, alpha = 0.1) +
      #geom_pointrange(aes(ymin = Lower, ymax = Upper)) +
      geom_line() +
      geom_hline(yintercept=1, linetype=2)+
      
      labs(x = "Fishing Year", y = "Index") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      scale_color_manual(values = cols)+
      scale_fill_manual(values = cols)+
      theme_cowplot() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        panel.grid.major = element_line(linetype = "dotted", colour = "grey", linewidth = 0.5)
      )
    
    
  } else {
    stan_labels <- c("TRUE" = "Standardised", "FALSE" = "Unstandardised")
    p <-  ggplot(index, aes(x = level, y = median, group = is_stan)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), 
                  data = filter(index, is_stan),
                  fill = fill, color = NA, alpha = 0.1) +
      geom_line(aes(linetype = is_stan, colour = is_stan)) +
      geom_point(aes(shape = is_stan, colour = is_stan), size = 3) +
      labs(x = "Fishing Year", y = "Index") +
      scale_shape_manual(values = c('TRUE' = 16, 'FALSE' = 1),
                         labels = stan_labels) +
      scale_colour_manual(values = c('TRUE' = fill, 'FALSE' = "grey40"),
                          labels = stan_labels) +
      scale_linetype_manual(
        values = c('TRUE' = "solid", 'FALSE' = "dashed"),
        labels = stan_labels) +
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
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
      ) +
      # For Binomial index we show both probabilities and relative index
      (if(component=='Binomial') {
        list( facet_wrap(~is_scaled, ncol = 1, scales = "free_y", 
                         labeller = as_labeller(c(
                           'FALSE' = "Binomial - Probabilities",
                           'TRUE'   = "Relative Index"
                         ))),
              theme(legend.position.inside = c(0.05, 0.35))
        )
        
        # For combined index we show all three indices
      } else if (component=="Combined"){
        # facet option  
        list(
          geom_hline(yintercept=1, linetype=2),
          facet_wrap(~Index, scales = 'free_y',ncol=1),
          scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))),
          guides(colour = "none", shape = "none", linetype = "none"),
          theme(panel.border = element_blank())
        )
        
      } else {
        NULL
      }) 
  }
  
  return(p)
}


