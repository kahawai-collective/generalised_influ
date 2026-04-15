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

get_DHARMAres <- function(fit) {
 
  # Extract original data
  raw_data <- fit$data

  if (is.null(raw_data) && inherits(fit, "survreg")) {
    raw_data <- eval(fit$call$data)
  }

  # Allow start_latitude and start_longitude in place of lat/lon
  if (!(("lat" %in% colnames(raw_data)) & ("lon" %in% colnames(raw_data)))) {
    raw_data$lat <- raw_data$start_latitude
    raw_data$lon <- raw_data$start_longitude
  }

  # Subset original data to predictors plus the coordinates even if they were not used as predictors
  raw_data <- subset(
    raw_data,
    select = union(all.vars(delete.response(terms(fit))), c('lat', 'lon'))
  )

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
    predictors = all.vars(delete.response(terms(fit)))
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

      # Create plotting dataframe
      plot_df <- data.frame(
        residuals = diag_metrics_df$scaledResiduals,
        predictor = pred_vec
      ) %>%
        add_count(predictor, name = "n_obs")

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
      }
      out$uniformity$p.value.cor = p.adjust(out$uniformity$p.value)
      # Low p-vals will have red boxplots
      box_colors <- ifelse(
        !is.na(out$uniformity$p.value.cor) & out$uniformity$p.value.cor < 0.05,
        "red",
        "grey30"
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
      ggplot(
        plot_df,
        aes(x = predictor, y = residuals, fill = predictor, alpha = n_obs)
      ) +
        geom_boxplot() +
        geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed") +
        scale_fill_manual(values = box_colors) +
        guides(fill = "none", alpha = "none") +
        labs(x = term_label, y = "DHARMa Residuals") +
        scale_x_discrete(drop = FALSE, labels = x_labels) +
        theme_bw() +
        dynamic_theme
    }),
    term_labels
  )
  return(plot_list)
}


#################################################################################################
# plot spatial residuals draft 15 APr
#################################################################################################

spatial_DHARMares <- function(diag_metrics, coastline = NULL, time_periods = 'all', grid_size=10, thresh=30) {

  # ---- prepare map data ------
data_sf <- st_as_sf(diag_metrics$diag_metrics, coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(3994) %>%
  st_make_valid()

grid_size <- grid_size*1000 # Define side length: convert km to meters for the EPSG:3994 projection

grid <- st_make_grid(data_sf, cellsize = grid_size, square = TRUE) %>%
  st_sf() %>%
  mutate(grid_id = row_number())

data_sf <- st_join(data_sf, grid, join = st_intersects)

bbox <- st_bbox(grid)
  
  if(!is.null(coastline)){
# Crop coastline to bbox
coastline <- coastline %>%
  st_transform(3994) %>%
  st_make_valid() %>%
  st_crop(bbox)
  }

  # ---- assign time periods to iterate through ------
if(identical(time_periods, "all")) time_periods <- sort(unique(data_sf$fyear))
  
  # ---- loop through time periods to generale plots ------

plots_list <- lapply(time_periods, function(period) {
  
  # Filter data to selected timeframe
  time_mask <- data_sf$fyear %in% period
  data_this_period <- data_sf[time_mask, ]
  
  # Recalculate dharma residulals filtering for time period and aggregating by grid cell
  dharma_recalculated <- recalculateResiduals(diag_metrics, sel = time_mask, group = data_sf$grid_id)
  
  # Get centroids for the spatial autocorrelation test
  # DHARMa needs coordinates for the groups
  grid_centers <- grid %>%
    filter(grid_id %in% data_this_period$grid_id) %>%
    st_centroid() %>%
    st_coordinates()  
  
  # Run Moran's I test
  Moran_test <- testSpatialAutocorrelation(
    dharma_recalculated, 
    x = grid_centers[,1], 
    y = grid_centers[,2],
    plot = F
  )
  p_val_text <- if(Moran_test$p.value < 0.001) "< 0.001" else round(Moran_test$p.value, 3)
  
  # bind together grid ids and DHARMa resids by grid
  res_df <- data.frame(
  grid_id = sort(unique(data_sf$grid_id[time_mask])),
  dharma_res = dharma_recalculated$scaledResiduals
)
    
# Join back to the grid for plotting
plot_grid <- grid %>%
  inner_join(res_df, by = "grid_id")

  
 if (length(period)>1) period_label <- paste(range(period), collapse = "-") else period_label <- period

# ----- plot code -----
 p <- ggplot() +
    geom_sf(data = plot_grid, aes(fill = dharma_res), color = NA) + 
   (if(!is.null(coastline)){
  geom_sf(data = coastline, fill = "grey80", color = "grey40") 
   }) +
    scale_fill_gradient2('DHARMa Residual',
                         low = "red",
                         mid = 'grey90',
                         high = "blue",
                         midpoint = 0.5, 
                        limits = c(0, 1)) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    theme_bw() +
    labs(title = paste0("Time period: ", period_label, 
                       "\n Moran's I p-val: ", p_val_text))
  
  
    # store p-value for further use e.g. in traffic light table    
p@meta$MoransIpval <- test_out$p.value
  return(p)

})
  
  wrap_plots(plots_list, ncol = 2) + 
  plot_layout(guides = 'collect') +
    plot_annotation(
    title = "DHARMa Moran's I test for distance-based autocorrelation",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  ) & 
  theme(legend.position = 'right')
  
}