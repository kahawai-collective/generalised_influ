#' Plot DHARMa Residuals for GLM2 and Survreg
#'
#' @description This function creates a standardized diagnostic plot for DHARMa residuals
#' using ggplot2.
#'
#' @param fit A fitted model object (of class glm, glm2, or survreg).
#' @param compare_fit A fitted model object if model needs to be compared with a previous model
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom stats dunif model.frame model.response
#' @import DHARMa
#' @importFrom cowplot plot_grid
#' @export
#'

get_DHARMAres <- function(fit, compare_fit  = NULL) {
 
  # Extract original data
  raw_data <- fit$data

  if (is.null(raw_data) && inherits(fit, "survreg")) {
    raw_data <- eval(fit$call$data)
  }

  

  # Identify lat and lon columns in a dataset (could be start_latitude, mean_latitude, etc)
  if (!(("lat" %in% colnames(raw_data)) & ("lon" %in% colnames(raw_data)))) {
    
    for (coord in c("lat", "lon")) {
            matches <- grep(coord, colnames(raw_data), ignore.case = TRUE, value = TRUE)
    
    if (length(matches) > 0) {
      raw_data[[coord]] <- raw_data[[matches[1]]]
      message("Notice: Auto-assigned '", matches[1], "' as the ", coord, " column.")
    } else {
      stop("Error: No ", coord, " column found in the dataset.")
    }
  }
  }

  # Subset original data to predictors plus the coordinates even if they were not used as predictors
  raw_data <- subset(
    raw_data,
    select = union(all.vars(delete.response(terms(fit))), c('lat', 'lon'))
  )

  # Identify new vessels and years in the data

  if(!is.null(compare_fit)) {
    compare_data <- compare_fit$data
    if (is.null(compare_data) && inherits(compare_fit, "survreg")) {
    compare_data <- eval(compare_fit$call$data)
  }
    # Identify which columns exist in both datasets
    common_cols <- intersect(names(raw_data), names(compare_data))
    
    new_factor_levels <- sapply(common_cols, function(col) {

      if(is.factor(raw_data[[col]])){
      setdiff(unique(raw_data[[col]]), unique(compare_data[[col]]))
      }
    }, simplify = FALSE)
    
    new_factor_levels <- Filter(function(x) length(x) > 0, new_factor_levels)
    }


  # Build dataframe to store diagnostics
  diag_metrics <- data.frame(
    # the reason for sub-setting to [1:nrow(data)] is to accommodate survreg which has extra column in response
    observed = as.vector(model.response(model.frame(fit)))[
      1:nrow(model.frame(fit))
    ],
    fitted = predict(fit, type = "response")
  )

  diag_metrics <- cbind(diag_metrics, raw_data)
  n <- nrow(diag_metrics)
  #_____________________________________________________________________________
  # Generate DHARMa Object
  #_____________________________________________________________________________

  if (inherits(fit, "survreg")) {
    simulatedResponse <- replicate(
      1000,
      rweibull(
        nrow(diag_metrics),
        shape = 1 / fit$scale,
        scale = predict(fit, type = "lp") %>% exp()
      )
    )

    DHARMa_obj <- DHARMa::createDHARMa(
      simulatedResponse = simulatedResponse,
      observedResponse = diag_metrics$observed,
      fittedPredictedResponse = diag_metrics$fitted,
      seed = 123
    )
  } else {
    DHARMa_obj <- DHARMa::simulateResiduals(
      fit,
      n = 1000,
      refit = FALSE,
      plot = FALSE,
      seed = 123
    )
  }

  # Extrac all recessary cells from DHARMa object

  diag_metrics$scaledResiduals <- residuals(DHARMa_obj)
  diag_metrics$observedResponse <- DHARMa_obj$observedResponse
  diag_metrics$fittedPredictedResponse <- DHARMa_obj$fittedPredictedResponse
  diag_metrics$simulatedResponse <- DHARMa_obj$simulatedResponse

  # Transform fitted values to ranks
  diag_metrics <- diag_metrics %>%
    mutate(
      rank_fitted = rank(fitted, ties.method = "average") / n,
      is_outlier = factor(
        scaledResiduals == 0 | scaledResiduals == 1,
        levels = c(FALSE, TRUE)
      )
    )
attr(diag_metrics, "predictors") <- all.vars(delete.response(terms(fit)))
  
  out <- list(
    scaledResiduals = diag_metrics$scaledResiduals,
    observedResponse = diag_metrics$observedResponse,
    fittedPredictedResponse = diag_metrics$fittedPredictedResponse,
    simulatedResponse = diag_metrics$simulatedResponse,
    refit = DHARMa_obj$refit,
    integerResponse = DHARMa_obj$integerResponse,
    diag_metrics = diag_metrics ,
    predictors = all.vars(delete.response(terms(fit))),
    new_factor_levels = if (exists("new_factor_levels")) new_factor_levels else NULL
  )

  class(out) <- "DHARMa"

  return(out)
}

#' Plot DHARMa Residuals for GLM2 and Survreg
#'
#' @description This function creates a standardized diagnostic plot for DHARMa residuals
#' using ggplot2.
#'
#' @param fit A fitted model object (of class glm, glm2, or survreg).
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom stats dunif model.frame model.response
#' @import DHARMa
#' @importFrom cowplot plot_grid
#' @export
#'
#'
plot_DHARMares <- function(diag_metrics) {

  diag_metrics <- diag_metrics$diag_metrics
  #______________________________________________________________________________
  #  Prepare residuals
  #______________________________________________________________________________

  # prepare quantile test and colour for quantile regression line
  q_test <- DHARMa::testQuantiles(
    diag_metrics$scaledResiduals,
    diag_metrics$rank_fitted,
    plot = F
  )
  p_vals <- q_test$pvals
  line_colours <- ifelse(p_vals <= 0.05, "red", "black")

  #______________________________________________________________________________
  #  Plotting code
  #______________________________________________________________________________

  # (1) ---Histogram of residuals versus uniform distribution---
  d1 <- ggplot(diag_metrics) +
    geom_histogram(
      aes(x = qnorm(scaledResiduals), y = after_stat(density)),
      stat = 'bin',
      binwidth = 0.05,
      fill = NA,
      colour = 'black'
    ) +
    labs(x = 'DHARMa residual', y = 'Density') +
    stat_function(
      fun = dnorm,
      #     args = list(min = 0, max = 1),
      color = "red",
      linewidth = 1
    ) +
    theme_classic() +
    theme(axis.title = element_text(size = 10))

  # (2) ---Residuals vs fitted values---

  # built-in function from DHARMa package, but result cannot be stored and used in patchwork
  # d2 <- DHARMa::plotResiduals(DHARMa, main = NULL)

  d2 <- ggplot(diag_metrics, aes(x = rank_fitted, y = scaledResiduals)) +
    geom_point(
      aes(shape = is_outlier, color = is_outlier, alpha = is_outlier),
      fill = 'white',
      stat = 'unique'
    ) +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 8)) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    scale_alpha_manual(values = c("FALSE" = 0.2, "TRUE" = 1)) +
    labs(x = 'Model predictions (rank transformed)', y = 'DHARMa residual') +
    geom_hline(
      yintercept = c(0.25, 0.5, 0.75),
      color = line_colours,
      linetype = "dashed",
      alpha = 0.5
    ) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, bs = "cs"),
      color = line_colours[2],
      se = FALSE
    ) +
    guides(shape = "none", color = "none", alpha = "none") +
    theme_classic() +
    theme(axis.title = element_text(size = 10))

  # (3) ---QQ plot---

  # built-in function from DHARMa package, but result cannot be stored and used in patchwork
  # d3 <- DHARMa::plotQQunif(DHARMa, testDispersion = FALSE,
  #                          testUniformity = FALSE,
  #                          testOutliers = FALSE, main = NULL)

  d3 <- ggplot(diag_metrics) +
    geom_qq(
      aes(sample = scaledResiduals),
      distribution = stats::qunif,
      dparams = list(min = 0, max = 1),
      cex = 0.3
    ) +
    labs(
      x = 'DHARMa resid. theoretical quantile',
      y = 'DHARMa resid. sample quantile'
    ) +
    geom_abline(intercept = 0, linetype = 'dotted', colour = 'blue') +
    theme_classic() +
    theme(axis.title = element_text(size = 10))

  # (4) ---Observed versus fitted values---
  d4 <- ggplot(diag_metrics) +
    geom_point(
      aes(fitted, observed),
      cex = 0.3,
      shape = 21,
      color = 'black',
      fill = 'white',
      stat = 'unique'
    ) +
    labs(x = 'Fitted value', y = 'Observed value') +
    geom_abline(intercept = 0, linetype = 'dotted', colour = 'blue') +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    theme_classic() +
    theme(axis.title = element_text(size = 10))

  p <- plot_grid(d1, d2, d3, d4, nrow = 2, align = 'hv')

  return(p)
}

#' Plot DHARMa Residuals Against Model Predictors
#'
#' @description
#' Plots DHARMa simulated residuals by levels of rpedictor variable.
#' Continuous variables are automatically binned. Boxplots are colored
#' red if the within-group Kolmogorov-Smirnov test for uniformity is significant
#' (p < 0.05, adjusted).
#'
#' @param fit A fitted model object (e.g., from \code{lme4}, \code{glmmTMB}, or \code{glm}).
#' @param model_data A \code{data.frame} containing the data used to fit the model.
#' @param simulationOutput An object of class \code{DHARMa}, typically generated by
#'   \code{\link[DHARMa]{simulateResiduals}}.
#'
#' @return A \code{patchwork} grid of \code{ggplot2} diagnostic plots.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_hline scale_fill_manual scale_alpha_continuous guides labs theme_bw theme element_text margin
#' @importFrom dplyr add_count
#' @importFrom magrittr %>%
#' @importFrom patchwork wrap_plots
#' @importFrom stats as.formula ks.test p.adjust na.omit
#'
#' @export
boxplot_DHARMares <- function(diag_metrics) {
  diag_metrics_df <- diag_metrics$diag_metrics
  term_labels <- diag_metrics$predictors

  plot_list <- setNames(
    lapply(term_labels, function(term_label) {
      # ----- Prepare data -----
      # Extract vector of predictor values
      pred_vec = diag_metrics_df[[term_label]]

      # Numeric terms are cut into factors
      if (
        is.numeric(pred_vec) &&
          (!all(pred_vec %% 1 == 0) | length(unique(pred_vec)) > 10)
      ) {
        breaks <- pretty(pred_vec, 30)
        step <- breaks[2] - breaks[1]
        breaks <- breaks[breaks <= max(pred_vec)]

        if (max(breaks) < max(pred_vec)) {
          breaks <- c(breaks, max(breaks) + step)
        }
        labels <- breaks[-length(breaks)] + (step / 2)
        pred_vec <- droplevels(cut(
          pred_vec,
          breaks = breaks,
          labels = labels,
          include.lowest = TRUE
        ))
      } else {
        pred_vec <- factor(pred_vec)
      }

       # Conduct kolmagorov-smirnov test.  Test taken directly from DHARMa package
      out = list()

      out$uniformity$details = suppressWarnings(by(
        diag_metrics_df$scaledResiduals,
        pred_vec,
        ks.test,
        'punif',
        simplify = TRUE
      ))
      out$uniformity$p.value = rep(NA, nlevels(pred_vec))
      for (i in 1:nlevels(pred_vec)) {
        out$uniformity$p.value[i] = out$uniformity$details[[i]]$p.value
        out$uniformity$KSstat[i] = out$uniformity$details[[i]]$statistic
      }
      out$uniformity$p.value.cor = p.adjust(out$uniformity$p.value)
      # Low p-vals will have red boxplots
      box_colors <- ifelse(
        !is.na(out$uniformity$p.value.cor) & out$uniformity$p.value.cor < 0.05,
        "red",
        "grey30"
      )

      # Create plotting dataframe
      plot_df <- data.frame(
        residuals = qnorm(diag_metrics_df$scaledResiduals),
        predictor = pred_vec
      ) %>%
        add_count(predictor, name = "n_obs") %>%
        mutate(
          is_new = as.character(predictor) %in% as.character(diag_metrics$new_factor_levels[[term_label]])
        )

      # ----- Plotting theme -----
      if ((length(levels(pred_vec)) > 12)) {
        # Vertical labels and tighter margins
        dynamic_theme <- theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.x.top = element_text(vjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 15))
        )
      } else {
        # Horizontal labels
        dynamic_theme <- theme(
          axis.text.x = element_text(angle = 0),
          axis.title.x = element_text(margin = margin(t = 10))
        )
      }

      # No x-labels for vessel_key
      x_labels <- if (grepl('vessel_key', term_label, ignore.case = TRUE)) {
        NULL
      } else {
        waiver()
      }

      # ----- Plot code -----
     p <-  ggplot(
        plot_df,
        aes(x = predictor, y = residuals, fill = predictor, alpha = n_obs)
      ) +
        geom_boxplot(outlier.shape = NA) +
        geom_hline(yintercept = c(-1.96, 0, 1.96), linetype = "dashed") +
        scale_fill_manual(values = box_colors) +
        guides(fill = "none", alpha = "none", linetype = 'none', linewidth = 'none') +
        labs(x = term_label, y = "DHARMa Residuals") +
        scale_x_discrete(drop = FALSE, labels = x_labels) +
        theme_bw() +
        dynamic_theme

      if (term_label %in% names(diag_metrics$new_factor_levels)) {
      p <- p + 
        aes(linewidth = is_new, linetype = is_new) +
        scale_linewidth_manual(values = c("TRUE" = 1.5, "FALSE" = 0.5)) +
        scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotdash"))
        
        # Calculate traffic light indicator:
        is_new       <- levels(pred_vec) %in% as.character(diag_metrics$new_factor_levels[[term_label]])
        meanKS_new   <- mean(out$uniformity$KSstat[is_new])
        meanKS_old   <- mean(out$uniformity$KSstat[!is_new])
        pct_increase <- (meanKS_new - meanKS_old) / meanKS_old 
        p@meta$indicatorKS  <- case_when(
          pct_increase <= 0.01 ~ '✔',
          pct_increase > 0.1 ~ '!',
          pct_increase <= 0.2 ~ '?'
        )
      }
      return(p)
    }),
    term_labels
  )
  return(plot_list)
}

#' Calculate and Plot Spatial DHARMa Residuals
#'
#' @description 
#' This function takes DHARMa diagnostic metrics, aggregates them onto a 
#' spatial grid, and calculates distance-based spatial autocorrelation (Moran's I) 
#' for specified time periods. It filters out grid cells with sparse data and outputs 
#' a combined `patchwork` spatial map of the residuals. 
#'
#' @param diag_metrics A list containing DHARMa diagnostic output (must include a nested 
#'   dataframe `$diag_metrics` with spatial coordinates `lon`, `lat` and time identifier `fyear`).
#' @param coastline An \code{sf} object representing the land/coastline to be masked over the plot. 
#'   Defaults to \code{NULL} (no coastline plotted).
#' @param plot_time_periods A vector, list of vectors, or the string \code{"all"}. Defines the time 
#'   blocks to iterate over. If \code{"all"} (default), iterates through all unique years.
#' @param grid_size Numeric. The side length of the spatial grid cells in kilometers. 
#'   Defaults to \code{10}.
#' @param thresh Numeric. The minimum number of observations required in a grid cell 
#'   for it to be included in the spatial analysis. Defaults to \code{3}.
#'
#' @return A \code{patchwork} plot object containing the combined spatial maps. A dataframe 
#'   summarising the Moran's I p-values for each time block is attached as an attribute 
#'   and can be accessed via \code{plot_object@meta$p_val_df}.
#' 
#' @import sf
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @importFrom DHARMa recalculateResiduals testSpatialAutocorrelation
#' @export
#'
#' 
spatial_DHARMares <- function(diag_metrics, coastline = NULL, plot_time_periods = 'all', grid_size=10, thresh=3, sea_only = FALSE) {

  # ---- prepare map data ------
  # Calculate extreme limits based on the data's distribution (use 3x IQR as a threshold)
  q_lat <- quantile(diag_metrics$diag_metrics$lat, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr_lat <- diff(q_lat)
  q_lon <- quantile(diag_metrics$diag_metrics$lon, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr_lon <- diff(q_lon)
  
  # Filter out the extreme spatial outliers
  data_sf <- diag_metrics$diag_metrics %>%
    mutate(orig_row_id = row_number()) %>%
    # drop_na(lat, lon) %>%
    filter(
      lat >= (q_lat[1] - 3 * iqr_lat) & lat <= (q_lat[2] + 3 * iqr_lat),
      lon >= (q_lon[1] - 3 * iqr_lon) & lon <= (q_lon[2] + 3 * iqr_lon)
    ) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
    st_transform(3994) %>%
    st_make_valid()


bbox <- st_bbox(data_sf)
  
  if(!is.null(coastline)){
# Crop coastline to bbox
coastline <- coastline %>%
  st_transform(3994) %>%
  st_make_valid() %>%
  st_crop(bbox)
     # Groom out data on land
    if(sea_only) {
      data_sf <- data_sf[lengths(st_intersects(data_sf, coastline)) == 0, ]
      }
  } else if (sea_only) {
    warning("sea_only = TRUE was requested, but no coastline object was provided. Skipping land filter.")
  }
  
  grid_size <- grid_size*1000 # Define side length: convert km to meters for the EPSG:3994 projection

grid <- st_make_grid(data_sf, cellsize = grid_size, square = TRUE) %>%
  st_sf() %>%
  mutate(grid_id = row_number())

 # Filtering for min recods here so that bbox zooms in on the data.  
data_sf <- st_join(data_sf, grid, join = st_intersects) %>%
  group_by(grid_id) %>%
  filter(n() >= thresh) %>%
  ungroup()
  
bbox <- st_bbox(data_sf)

  # ---- Assign time periods to iterate through ------
if(identical(plot_time_periods, "all")) plot_time_periods <- sort(unique(data_sf$fyear))
  years_for_stats <- as.list(tail(sort(unique(as.numeric(as.character(data_sf$fyear)))), 10))

  time_periods <- union(years_for_stats, plot_time_periods)
  
  # ---- Loop through time periods to generale plots ------

results_list <- lapply(time_periods, function(period) {
  
  # Filter data to selected timeframe and metting the threshold

  data_this_period <- data_sf %>%
    filter(fyear %in% period) %>%
    group_by(grid_id) %>%
    filter(n() >= thresh) %>%
    ungroup()

  # Label formatting
 period_label <- if (length(period)>1) paste(range(period), collapse = "-") else period
  
    if(nrow(data_this_period) == 0) {
      warning(paste("Time period", period_label, "dropped: No grid cells met the threshold of", thresh))
      return(NULL) 
    }
  
  # Recalculate DHARMa residulals, while filtering for time period, minimum obs, and aggregating by grid cell
  # This function needs length(group) to be the same as length of simulated residuals in diag_metrics
  # Create an empty vector the exact length of the original DHARMa data
  full_length_group <- rep(NA, length(diag_metrics$scaledResiduals))
  
  # 2. Inject the valid grid IDs into their exact original row positions
  full_length_group[data_this_period$orig_row_id] <- data_this_period$grid_id

  dharma_recalculated <- recalculateResiduals(diag_metrics, sel = data_this_period$orig_row_id,   
    group = full_length_group)
  
    
  # Run Moran's I test
    # Get centroids for the spatial autocorrelation test
    # DHARMa needs coordinates for the groups

  grid_centers <- grid %>%
    filter(grid_id %in% data_this_period$grid_id) %>%
    st_centroid() %>%
    st_coordinates()  

  Moran_test <- testSpatialAutocorrelation(
    dharma_recalculated, 
    x = grid_centers[,1], 
    y = grid_centers[,2],
    plot = F
  )

  p_val_text <- if(Moran_test$p.value < 0.001) "< 0.001" else round(Moran_test$p.value, 3)

  # bind together grid ids and DHARMa resids by grid
  res_df <- data.frame(
  grid_id = sort(unique(data_this_period$grid_id)),
  dharma_norm = scales::squish(qnorm(dharma_recalculated$scaledResiduals), c(-4, 4), only.finite = FALSE)
)
  

  p <- NULL
  stats_row <- NULL
  
  # ---------------------------------------------
  # OPTION A: Generate Plot (if in plot_time_periods)
  # ---------------------------------------------

    if ( list(period) %in% plot_time_periods){
# Join back to the grid for plotting
plot_grid <- grid %>%
  inner_join(res_df, by = "grid_id")

# ----- plot code -----
 p <- ggplot() +
    geom_sf(data = plot_grid, 
      aes(fill = dharma_norm),
      color = NA) + 
   (if(!is.null(coastline)){
  geom_sf(data = coastline, fill = "grey80", color = "grey40") 
   }) +
    scale_fill_gradientn(
    name = 'Residual',
    colors = c("darkred", "tomato", "grey90","cornflowerblue", "darkblue"),
    values = scales::rescale(c(-4, -2, 0, 2, 4)),
    limits = c(-4, 4)    
  ) +
   scale_x_continuous(n.breaks = 4) +
   scale_y_continuous(n.breaks = 4) +
   coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
   labs(title = paste0("Time period: ", period_label, 
                       "\nMoran's I p-val: ", p_val_text)) +
   theme_bw() +
   theme(
      plot.title = element_text(size = 8, lineheight = 1.1, margin = margin(b = 5)),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8), 
      axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8) 
    )
  }

  # ---------------------------------------------
  # OPTION B: Generate Stats (if in years_for_stats)
  # ---------------------------------------------

  if (list(period) %in% years_for_stats) {
  # store Moran's I p-value for further use e.g. in traffic light table 
  # Rescale Moran's I to Z value so that we can compare it across years
  MoransI_z_score  <- as.numeric((Moran_test$statistic["observed"] - Moran_test$statistic["expected"]) / Moran_test$statistic["sd"])

  # calculate proportion of spatial residuals in the tails of the distribution
  tail_prop <- mean(abs(res_df$dharma_norm) > 1.96)
    
    stats_row <- data.frame("period" = period_label,
    'MoransI_z_score' = MoransI_z_score,
    'tail_prop' = tail_prop
    )
  }
 
  return(list(plot = p, stats_row = stats_row))
}
)
  plots_list <- lapply(results_list, function(x) x$plot)
  plots_list <- Filter(Negate(is.null), plots_list)

  # ---- Combine plots ------

  combined_plot <- wrap_plots(plots_list, ncol = 2) + 
  plot_layout(guides = 'collect')  & 
  theme(legend.position = 'right')
    
  # ---- Extract Metadata Table ------
 stats_df <- bind_rows(lapply(results_list, function(x) x$stats_row))
  
  combined_plot@meta$stats <- stats_df

  return (combined_plot)
}

