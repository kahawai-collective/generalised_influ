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
#' @importFrom cowplot theme_cowplot 
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



#' Plot multiple CPUE indices for comarison
#'
#' Standardises and visualises multiple CPUE time series on a single plot. Supports 
#' geometric mean normalisation over a common overlap period or a -1 to 1 
#' min-max scaling for ENSO-style comparisons.
#'
#' @param cidx A named list of data frames containing the primary CPUE indices.
#' @param CPUE_set Character vector. The names of the elements in \code{cidx} to be plotted.
#' @param component Character vector. The specific \code{Index} type (e.g., 'Positive', 
#'   'Binomial', 'Combined') to extract from each set. Defaults to 'Positive'.
#' @param alt_CPUE1 Optional data frame. An additional index to include in the comparison.
#' @param componentAlt1 Character string. The \code{Index} type to extract from \code{alt_CPUE1}.
#' @param series_alt Character string. The legend label for \code{alt_CPUE1}.
#' @param alt_CPUE2 Optional data frame. A second additional index to include.
#' @param componentAlt2 Character string. The \code{Index} type to extract from \code{alt_CPUE2}.
#' @param series_alt2 Character string. The legend label for \code{alt_CPUE2}.
#' @param normalise_ENSO Logical. If \code{TRUE}, uses min-max scaling to bound data 
#'   between -1 and 1. If \code{FALSE}, uses geometric mean scaling relative to the 
#'   overlap period. Defaults to \code{FALSE}.
#' @param uncert Logical. If \code{TRUE}, adds uncertainty bars.
#'
#' @details 
#' The function identifies the "overlap period"—the years common to all included 
#' series—and calculates a geometric mean for each series based only on those years. 
#' This ensures that the relative trends are comparable even if the time series 
#' have different start or end years.
#' 
#' @return A \code{ggplot2} object.
#' 
#' @import dplyr
#' @import ggplot2
#' @importFrom cowplot theme_cowplot 
#' @export

compare_indices <- function(cidx, 
                            CPUE_set, 
                            component = 'Positive', 
                            alt_CPUE1 = NULL, 
                            componentAlt1=component[1],
                            series_alt1 = NULL,
                            alt_CPUE2 = NULL, 
                            componentAlt2=component[1],
                            series_alt2 = NULL,
                            normalise_ENSO = FALSE,
                            uncert=F){
  
  if(length(component)==1) component <- rep(component, length(CPUE_set))
  
  # helper function to filter idx and add idx rescaled between -1 and 1.
  process_idx <- function(idx, target_index) {
    idx %>% 
      filter(is_stan, is_scaled, Index %in% target_index) %>%
      mutate(level = as.numeric(as.character(level)),
             index.norm = 2 * ((median - min(median)) / (max(median) - min(median))) - 1)
  }
  
  indices <- setNames(lapply(seq_along(CPUE_set), function(l) {
    process_idx(cidx[[CPUE_set[l]]], component[l])
  }), names(cidx[CPUE_set]))
  
  
  joint <- bind_rows(indices, .id = 'Series')
  
  
  if(!is.null(alt_CPUE1)) joint <- bind_rows(joint, process_idx(alt_CPUE1, componentAlt1) %>% mutate(Series = series_alt1))
  if(!is.null(alt_CPUE2)) joint <- bind_rows(joint, process_idx(alt_CPUE2, componentAlt2) %>% mutate(Series = series_alt2))
  
  
  overlap <- joint %>%
    group_by(Series, Index) %>%
    summarise(min_y = min(level),
              max_y = max(level), 
              .groups = "drop") %>%
    summarise(start = max(min_y), 
              end = min(max_y))
  
  joint <- joint %>%
    group_by(Series, Index) %>%
    mutate(
      # Calculate gmean only on the subset that falls in the overlap
      gmeans = gmean(median[level >= overlap$start & level <= overlap$end]),
      index  = median / gmeans,
      Lower = Lower / gmeans,
      Upper = Upper / gmeans
    ) %>%
    rename(`Index type` = Index)
  
  
  y_col_name <- if (normalise_ENSO)  "index.norm" else "index"
  myColors <- c('dodgerblue', '#F5B915FF', '#08235FFF', '#4D9221',  "purple4" ,  "violetred")
  n_series <- length(unique(joint$Series))
  
  
  g <- ggplot(joint, aes(level, 
                         y = .data[[y_col_name]],
                         ymin = Lower,
                         ymax = Upper,
                         linetype = `Index type`,
                         shape = `Index type`,
                         col = Series, 
                         group = interaction(Series, `Index type`))) +
    
    geom_line() +
    geom_point() +
    # Conditional Layers
    ( if (normalise_ENSO){
      scale_y_continuous("CPUE index")
    } else {
      list(
        scale_y_continuous("CPUE index", limits = c(0, NA)),
        geom_hline(yintercept=1, linetype=3),
        if(uncert)  geom_linerange()
      )
    }) +
    
    scale_x_continuous("Fishing year", breaks = unique(joint$level)) +
    scale_color_manual(values = myColors) +
    theme_cowplot() +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45),
          panel.grid.major = element_line(colour = "grey90"),
          legend.position = "bottom",
          legend.direction = "vertical")
  
  
  if (length(unique(joint$`Index type`)) == 1) {
    # number of cols is 2 for 4+ series, and one row otherwise
    g <- g +guides(linetype = "none", shape = "none", colour = guide_legend(ncol = min(n_series, 2 + (n_series == 3))))
    
  }
  
  return(g)
  
}
