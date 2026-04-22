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
                            # component = 'Positive', 
                            alt_CPUE1 = NULL, 
                            # componentAlt1=component[1],
                            series_alt1 = NULL,
                            alt_CPUE2 = NULL, 
                            # componentAlt2=component[1],
                            series_alt2 = NULL,
                            normalise_ENSO = FALSE,
                            uncert=F){
  
  # if(length(component)==1) component <- rep(component, length(CPUE_set))
  
  # helper function to filter idx and add idx rescaled between -1 and 1.
process_idx <- function(idx, series_set) {
    idx %>% 
      filter(is_stan, is_scaled, Series %in% series_set, tolower(Index)==selected_idx) %>%
      group_by(Series) %>%
      mutate(index.norm = 2 * ((median - min(median)) / (max(median) - min(median))) - 1)
  }
    
  indices <- process_idx(cidx, CPUE_set)
  
  if(!is.null(alt_CPUE1)) indices <- bind_rows(indices %>%
    mutate(Series = paste(Series, 'UPDATE')), process_idx(alt_CPUE1, series_alt1))%>%
    arrange(Series)
    
  
  if(!is.null(alt_CPUE2)) indices <- bind_rows(indices %>% 
    mutate(Series = paste(Series, 'UPDATE_2')), process_idx(alt_CPUE2, series_alt2) )%>%
    arrange(Series)
    
  overlap <- indices %>%
    group_by(Series, Index) %>%
    summarise(min_y = min(level),
              max_y = max(level), 
              .groups = "drop") %>%
    summarise(start = max(min_y), 
              end = min(max_y))
  
  indices <- indices %>%
    group_by(Series, Index) %>%
    mutate(
      # Calculate gmean only on the subset that falls in the overlap
      gmeans = gmean(median[level >= overlap$start & level <= overlap$end]),
      index  = median / gmeans,
      Lower = Lower / gmeans,
      Upper = Upper / gmeans
    ) %>%
    rename(`Index type` = Index)
  
  
  #  ----- Time series stats -----
  # _______________________________________________________________________________
if(!is.null(alt_CPUE1)) {
  

level_divergence = function(current, last) {
  1-sum(abs(current-last), na.rm = T)/sum(pmax(current, last), na.rm = T)
}

trend_divergence <- function(current, last, level, mode = "overlap") {
  
  df <- data.frame(cur = current, lst = last, lvl = level)
  
  # If calculating for overlapping period: Removes any row where either series is NA 
  if (mode == "overlap")   df <- na.omit(df)
     
  smooth_cur  <- predict(loess(cur ~ lvl, data = df, span = 0.9, control = loess.control(surface = "direct"), na.action = na.omit))
  smooth_last <- predict(loess(lst ~ lvl, data = df, span = 0.9, control = loess.control(surface = "direct"), na.action = na.omit),
   newdata = data.frame(lvl = df$lvl))
  
  d_cur  <- diff(smooth_cur)
  d_last <- diff(smooth_last)
  
  # Create exponential weights. We want the last year to be the heaviest.
  weights <- exp(seq(0, 1, length.out = length(d_cur)))

  # Compare the "energy" and direction of the slopes
  1 - sum(weights*abs(d_cur - d_last), na.rm = TRUE) / (sum(weights*pmax(abs(d_cur), abs(d_last)), na.rm = TRUE)+0.01)

  }

  series_stats <- indices %>%
    ungroup() %>%
    select(level, Series, index) %>%
          
    pivot_wider(names_from = Series, values_from = index) %>%
    rename_with(~ ifelse(grepl("UPDATE", .x), "current", "last"), 
                .cols = -level) %>%
    arrange(level) %>%
    summarise(
      Level_Div = level_divergence(current, last),
      Trend_Div_Overlap = trend_divergence(current, last, level, mode = "overlap"),
      Trend_Div_Full = trend_divergence(current, last, level, mode = "full")
    )

} else {series_stats <- NULL}
  
 #__________________________________________________________________________________________________________
  # End of series stats code 
  
  y_col_name <- if (normalise_ENSO)  "index.norm" else "index"
  myColors <- c( '#F5B915FF', 'dodgerblue', '#08235FFF', '#4D9221',  "purple4" ,  "violetred")
  n_series <- length(unique(indices$Series))
  
  
  g <- ggplot(indices, aes(level, 
                           y = .data[[y_col_name]],
                           ymin = Lower,
                           ymax = Upper,
                           linetype = `Index type`,
                           shape = `Index type`,
                           col = Series, 
                           group = interaction(Series, `Index type`))) +
    # Conditional Layers
    ( if (normalise_ENSO){
      scale_y_continuous("CPUE index", limits = c(-1.1, 1.1))
    } else {
      list(
        scale_y_continuous("CPUE index", limits = c(0, NA)),
        geom_hline(yintercept=1, linetype=3),
        if(uncert)  geom_linerange()
      )
    }) +
    geom_line() +
    geom_point() +
    scale_x_continuous("Fishing year", breaks = unique(indices$level)) +
    scale_color_manual(values = myColors) +
    theme_cowplot() +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45),
          panel.grid.major = element_line(colour = "grey90"),
          legend.position = "bottom",
          legend.direction = "vertical")
  
  
  if (length(unique(indices$`Index type`)) == 1) {
    # number of cols is 2 for 4+ series, and one row otherwise
    g <- g +guides(linetype = "none", shape = "none", colour = guide_legend(ncol = min(n_series, 2 + (n_series == 3))))
    
  }
  
  g@meta$series_stats <- series_stats

  return(g)
  
  }


#' Plot Status of the stock 
#'
#' Generates a multi-panel plot showing CPUE indices relative to reference levels (if exist), 
#' removals, and relative exploitation rates.
#'
#' @param cidx A list of CPUE index data frames generated with get_index().
#' @param CPUE_set Character vector of series names to include.
#' @param component Character vector of index components, one of c('BInomial', 'Positive', 'Combined').
#' @param ref_period Numeric vector of years for the reference period.
#' @param period_type Reference type one of c('target', 'soft_limit', 'hard_limit').
#' @param ref_series Index of the reference series in `CPUE_set`.
#' @param ref_index The name of the reference index component, one of c('BInomial', 'Positive', 'Combined').
#' @param landings_data Data frame of landings/removals.
#' @param plot_exploitation Logical; if TRUE, includes the exploitation rate panel.
#' @param cpue_smooth Logical; applies loess smoothing to the CPUE.
#' @param bmsy_proxy Numeric proxy value (default 40).
#'
#' @return A patchwork ggplot object.
#' @importFrom cowplot theme_cowplot 
#' @export

plot_sos <- function(cidx, 
                     CPUE_set, 
                     ref_period=NULL, 
                     period_type='target', 
                     ref_series = 1, 
                     # ref_index = 'Positive', 
                     landings_data = NULL, 
                     plot_exploitation = TRUE, 
                     cpue_smooth = FALSE, 
                     bmsy_proxy = 40){
  
  
  ref_name <- unique(cidx$Series)[ref_series]
  
  indices <- cidx %>% 
    filter(is_stan, is_scaled, Series %in% CPUE_set, tolower(Index)==selected_idx) %>%
    mutate(is_ref = (Series == ref_name))
  
  overlap <- indices %>%
    group_by(Series, Index) %>%
    summarise(min_y = min(level),
              max_y = max(level), 
              .groups = "drop") %>%
    summarise(start = max(min_y), 
              end = min(max_y))
  
  indices <- indices %>%
    group_by(Series, Index) %>%
    mutate(
      # Calculate gmean only on the subset that falls in the overlap
      gmeans = gmean(median[level >= overlap$start & level <= overlap$end]),
      index  = median / gmeans,
      Lower  = Lower / gmeans,
      Upper  = Upper / gmeans
    ) 
  #_____________________________________________________________________________
  # Stock status plot
  #_____________________________________________________________________________
  
  g1 <- ggplot(indices, aes(level, index, group = interaction(Series,Index), linetype = Index, color = is_ref))
  
  
  if(!is.null(ref_period)){
    
    # Extract reference index
    this_idx <-indices[indices$is_ref,]
    
    # Smooth it if desired, subset to ref period only
    if (cpue_smooth) {
      ll      <- loess(index~level, data=this_idx)
      ref_vals <- ll$fitted[this_idx$level %in% ref_period]
    } else {
      ref_vals <- this_idx$index[this_idx$level %in% ref_period]
    }
    
    ref_mult <- case_when(
      period_type == 'target' ~ 1,
      period_type == 'soft_limit' ~ bmsy_proxy/20,
      period_type == 'hard_limit' ~ bmsy_proxy/10
    )
    
    # Calculate base value
    b <- gmean(ref_vals) * ref_mult
    
    # Add B10, B20, B40 and reference period to the plot
    g1 <- g1 + geom_hline(yintercept = b * c(1, 20/bmsy_proxy, 10/bmsy_proxy), 
                          linetype = c(5,2,4), col = c("seagreen", "orange", "tomato")) +
      geom_vline(xintercept = range(ref_period), linetype = 'dashed', col = 'blue')
    
  }
  
  # Add indices to the plot  
  g1 <- g1 +
    
    # Layer 1: Background (All non-reference series)
    geom_line(data = filter(indices, !is_ref), color = "grey80") +
    geom_point()+
    
    # Layer 2: Foreground (The reference series only)
    geom_line(data = filter(indices, is_ref), color = "black") +
    geom_linerange(data = filter(indices, is_ref), aes(ymin = Lower, ymax = Upper), size = 0.5) +
    
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey80"), guide = "none") +
    scale_y_continuous('CPUE index', limits = c(0, NA)) +
    theme_cowplot() +
    theme(panel.grid.major = element_line(colour = 'grey90')) +
    (if(plot_exploitation){
      list(
        xlab('') ,   
        theme(axis.text.x = element_blank())
      )
    } else{
      list(
        xlab('Fishing year') , 
        theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = rel(0.7)))
      )
    })
  
  if (length(unique(indices$Index)) == 1) {g1 <- g1 +  scale_linetype(guide = 'none') }
  
  if (cpue_smooth) {g1 <- g1 + geom_smooth(data = this_idx, aes(x=level, y=index))}
  
  # If no landings data, we exit here
  if(is.null(landings_data)) return(g1)  
  
  #_____________________________________________________________________________
  # If have landings data: Removals plot and Exploitation rate plot
  #_____________________________________________________________________________
  
  landings <- landings_data %>% filter(destination_type != 'X') %>%
    group_by(level = as.integer(fishyear)) %>%
    summarise(landings = sum(gwt, na.rm=T)/1000) %>%
    inner_join(indices, by = c('level')) %>%
    mutate(erate = landings/index,
           erate = erate/gmean(erate))
  
  # --- Removals plot ---
  
  g2 <- ggplot(landings,aes(level, landings))  +
    scale_y_continuous('Removals (t)', limits=c(0, NA)) +
    geom_line()+
    geom_point()+
    xlab('') +
    theme_cowplot() +
    theme(axis.text.x = element_blank())
  
  # --- Exploitation rate plot ---
  if (plot_exploitation) {
    
    
    g3 <- ggplot(data = landings %>% filter(is_ref), aes(x=level, y=erate))  +
      # scale_x_continuous('Fishing year',breaks=unique(joint$level),limits=range(unique(joint$level))) +
      scale_x_continuous('Fishing year')+
      scale_y_continuous('Relative exploitation rate', limits = c(0,NA)) +
      geom_line()+
      geom_point()+
      theme_cowplot() +
      theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = rel(0.7))) +
      if(!is.null(ref_period)){
        geom_hline(yintercept=gmean(landings$erate[ landings$is_ref & landings$level %in% ref_period])/ref_mult,
                   linetype='longdash', col='seagreen')
      }
    
    # three plots
    plot_combined <- (g2 / g1 / g3) + 
      plot_layout(guides = "collect") + 
      plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) & 
      theme(panel.grid.major = element_line(colour = 'grey90'),
            plot.tag.position = c(0.06, 1.15),
            plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
    
    # two plots
  } else {
    plot_combined <- (g2 / g1) + 
      plot_layout(guides = "collect") + 
      plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) & 
      theme(panel.grid.major = element_line(colour = 'grey90'),
            plot.tag.position = c(0.06, 1.1),
            plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
  }
  
  
  return(plot_combined)
}
