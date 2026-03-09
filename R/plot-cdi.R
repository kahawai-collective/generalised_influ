#' Bayesian version of the CDI plot
#' 
#' The CDI plot presents the coefficients for the variable of interest (top-left panel), the spread of the data 
#' (bottom-left panel), and the influence statistic (bottom-right panel).
#' 
#' @param fit An object of class \code{brmsfit}.
#' @param xfocus The column name of the variable to be plotted on the x axis. 
#'   This column name must match one of the column names in the 
#'   \code{data.frame} that was passed to \code{brm} as the \code{data} argument.
#' @param yfocus The column name of the variable to be plotted on the y axis. 
#'   This column name must match one of the column names in the 
#'   \code{data.frame} that was passed to \code{brm} as the \code{data} argument. 
#'   This is generally the temporal variable in a generalised linear model (e.g. year).
#' @param hurdle If a hurdle model then use the hurdle.
#' @param sort_coefs Should the coefficients be sorted from highest to lowest.
#' @param axis.text.x.bl Include the x axis labels on the bottom-left (bl) bubble plot panel.
#' @param xlab The x axis label.
#' @param ylab The y axis label.
#' @param colour The colour to use in the plot.
#' @param p_margin The margin between panels on the plot. This is passed to \code{margin} within \code{theme}.
#' @param legend To show the legend or not.
#' @param sum_by Sum to 1 by row, sum to 1 by column, sum to 1 across all data, or raw. The size of the bubbles will be 
#'   the same for all and raw, but the legend will change from numbers of records to a proportion.
#' @param ... Further arguments passed to nothing.
#' @return a \code{ggplot} object.
#' @seealso \code{\link{get_coefs}}, \code{\link{get_influ}}, \code{\link{plot_bubble}}
#' @importFrom gtable is.gtable gtable_filter
#' @importFrom stats poly
#' @importFrom tidyselect all_of
#' @importFrom stats median
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @export
#' 
plot_bayesian_cdi <- function(fit, 
                              xfocus = "area", yfocus = "fishing_year",
                              xlab = NULL, ylab = NULL, 
                              hurdle = FALSE,
                              sort_coefs = FALSE, 
                              axis.text.x.bl = TRUE,
                              colour = "purple", 
                              p_margin = 0.05, 
                              legend = TRUE, 
                              sum_by = "row", ...) {
  
  if (!is.brmsfit(fit)) stop("fit is not an object of class brmsfit.")
  if (is.null(xlab)) xlab <- xfocus
  if (is.null(ylab)) ylab <- yfocus
  y_coefs <- "Conditional effect"
  
  # Identify the type of variable we are dealing with
  type <- id_var_type(fit = fit, xfocus = xfocus, hurdle = hurdle)
  
  # Posterior samples of coefficients
  # if (type %in% c("random_effect")) {
  if (type %in% c("fixed_effect", "random_effect")) {
    coefs <- get_coefs(fit = fit, var = xfocus, hurdle = hurdle)
  } else {
    # this would plot the marginal/conditional effect, but if it is a hurdle model it ignores the hurdle bit
    coefs <- get_marginal(fit = fit, var = xfocus)# %>%
    # mutate(value = log(value))
    # library(marginaleffects)
    # pred <- predictions(model = fit)
    # coefs <- get_coefs_raw(fit = fit, var = xfocus)
  }
  
  # If using the lognormal distribution then transform the coefs
  if (fit$family$family == "lognormal") {
    coefs <- coefs %>% mutate(value = exp(.data$value))
    y_coefs <- "Coefficient"
  }
  
  # Model data
  if (is.numeric(coefs$variable)) {
    data <- fit$data %>% select(all_of(c(yfocus, xfocus)))
    length.out <- 15
    dmin <- min(data[,xfocus])
    dmax <- max(data[,xfocus])
    breaks <- seq(dmin, dmax, length.out = length.out)
    midpoints <- breaks[-length(breaks)] + diff(breaks) / 2
    data[,xfocus] <- cut(data[,xfocus], breaks = breaks, labels = sprintf("%.2f", round(midpoints, 2)), include.lowest = TRUE)
  } else {
    data <- fit$data
  }
  
  # Sort the coefficients if required
  sort_order <- NULL
  if (sort_coefs) {
    coefs_1 <- coefs %>%
      group_by(.data$variable) %>%
      summarise(value = median(.data$value))
    coefs_s <- coefs_1 %>%
      arrange(.data$value) %>% 
      select(.data$variable)
    
    # reorder coefficients
    coefs$variable <- factor(coefs$variable, levels = coefs_s$variable)
    
    # reorder bubbles
    coefs_o <- match(coefs_1$variable, coefs_s$variable)
    bubble_o <- data.frame(xfocus = data[,xfocus], order = as.numeric(data[,xfocus])) %>%
      distinct()
    bubble_o <- bubble_o[match(coefs_o, bubble_o$order),]
    data[,xfocus] <- factor(data[,xfocus], levels = bubble_o$xfocus)
    sort_order <- bubble_o$xfocus
  }
  
  # Influence
  influ <- get_influ2(fit = fit, group = c(yfocus, xfocus), hurdle = hurdle)
  
  # Extract the legend on its own
  g2 <- function(a.gplot) {
    if (!is.gtable(a.gplot))
      a.gplot <- ggplotGrob(a.gplot)
    gtable_filter(a.gplot, 'guide-box', fixed = TRUE)
  }
  
  # The bubble plot (bottom-left) and the legend for the bubble plot (top-right)
  p3a <- plot_bubble(df = data, group = c(yfocus, xfocus), sum_by = sum_by, 
                     xlab = xlab, ylab = ylab, zlab = "", fill = colour, sort_order = sort_order)
  
  p2 <- g2(p3a)
  
  if (axis.text.x.bl) {
    p3 <- p3a + theme(legend.position = "none", plot.margin = margin(t = p_margin, r = p_margin, unit = "cm"), 
                      axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    p3 <- p3a + theme(legend.position = "none", plot.margin = margin(t = p_margin, r = p_margin, unit = "cm"), 
                      axis.text.x = element_blank())
  }
  
  # The coefficients (top-left)
  p1 <- ggplot(data = coefs, aes(x = .data$variable, y = .data$value)) +
    labs(x = NULL, y = y_coefs) +
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(b = p_margin, r = p_margin, unit = "cm"))
  
  if (is.numeric(coefs$variable)) {
    p3 <- p3 + scale_x_discrete(expand = expansion(mult = 0.05))
    
    p1 <- p1 +
      stat_summary(geom = "ribbon", alpha = 0.5, fill = colour, 
                   fun.min = function(x) quantile(x, probs = 0.025), 
                   fun.max = function(x) quantile(x, probs = 0.975)) +
      stat_summary(fun = "median", geom = "line", colour = colour) +
      scale_x_continuous(position = "top", breaks = midpoints, minor_breaks = NULL, expand = expansion(mult = 0.05)) +
      coord_cartesian(xlim = c(midpoints[1], midpoints[length(midpoints)]))
  } else {
    p1 <- p1 +
      # geom_point() +
      geom_violin(colour = colour, fill = colour, alpha = 0.5, draw_quantiles = 0.5, scale = "width") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_x_discrete(position = "top")# +
    # scale_x_discrete(position = "top", breaks = midpoints, minor_breaks = NULL, expand = expansion(mult = 0.05)) +
    # coord_cartesian(xlim = c(midpoints[1], midpoints[length(midpoints)]))
  }
  
  # The influence plot (bottom-right)
  p4 <- ggplot(data = influ, aes_string(x = as.character(yfocus))) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_violin(aes(y = exp(.data$delta)), colour = colour, fill = colour, alpha = 0.5, draw_quantiles = 0.5, scale = "width") +
    # geom_violin(aes(y = .data$delta), colour = colour, fill = colour, alpha = 0.5, draw_quantiles = 0.5, scale = "width") +
    coord_flip() +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = "Influence") +
    theme_bw() +
    theme(legend.position = "none", plot.margin = margin(t = p_margin, l = p_margin, unit = "cm"))
  
  if (legend) {
    p <- p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2, heights = c(1, 2), widths = c(2, 1))
  } else {
    pv <- ggplot() + theme_void()
    p <- p1 + pv + p3 + p4 + plot_layout(nrow = 2, ncol = 2, heights = c(1, 2), widths = c(2, 1))
  }
  
  return(p)
}


#' Bayesian version of the CDI plot (depreciated)
#' 
#' @param fit a model fit
#' @param xfocus The column name of the variable to be plotted on the x axis. This column name must match one of the
#'   column names in the \code{data.frame} that was passed to \code{brm} as the \code{data} argument.
#' @param yfocus The column name of the variable to be plotted on the y axis. This column name must match one of the
#'   column names in the \code{data.frame} that was passed to \code{brm} as the \code{data} argument. This is generally the
#'   temporal variable in a generalised linear model (e.g. year).
#' @param hurdle if a hurdle model then use the hurdle
#' @param xlab the x axis label
#' @param ylab the y axis label
#' @param colour the colour to use in the plot
#' @return a ggplot object
#' 
#' @importFrom gtable is.gtable gtable_filter
#' @importFrom stats poly
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @export
#' 
plot_bayesian_cdi2 <- function(fit,
                               xfocus = "area", yfocus = "fishing_year",
                               hurdle = FALSE,
                               xlab = "Month", 
                               ylab = "Fishing year", 
                               colour = "purple") {
  
  # Posterior samples of coefficients
  coefs <- get_coefs(fit = fit, var = xfocus, normalise = TRUE, hurdle = hurdle)
  n_iterations <- max(coefs$iteration)
  
  get_midpoint <- function(cut_label) {
    mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ","))))
  }
  
  # Model data
  is_poly <- FALSE
  if (any(grepl("poly", coefs$variable))) {
    is_poly <- TRUE
    data <- fit$data %>%
      select(-starts_with("poly"))
    dmin <- min(data[,xfocus])
    dmax <- max(data[,xfocus])
    data[,xfocus] <- cut(data[,xfocus], breaks = seq(dmin, dmax, length.out = 20), include.lowest = TRUE)
    # breaks <- unique(quantile(data[,xfocus], probs = seq(0, 1, length.out = 15)))
    # data[,xfocus] <- cut(data[,xfocus], breaks = breaks, include.lowest = TRUE)
    data[,xfocus] <- sapply(data[,xfocus], get_midpoint)
    
    z <- poly(fit$data[,xfocus], 3)
    x_new <- data.frame(id = 1:length(unique(data[,xfocus])), variable = sort(unique(data[,xfocus])))
    x_poly <- poly(x_new$variable, 3, coefs = attr(z, "coefs"))
    
    # Do the matrix multiplication
    Xbeta <- matrix(NA, nrow = n_iterations, ncol = nrow(x_poly))
    for (i in 1:n_iterations) {
      Xbeta[i,] <- x_poly %*% filter(coefs, .data$iteration == i)$value
    }
    coefs <- melt(Xbeta, varnames = c("iteration", "id")) %>%
      left_join(x_new, by = "id") %>%
      select(-id)
  } else if (length(unique(coefs$variable)) == 1) {
    data <- fit$data# %>%
    # select(xfocus)
    dmin <- min(data[,xfocus])
    dmax <- max(data[,xfocus])
    data[,xfocus] <- cut(data[,xfocus], breaks = seq(dmin, dmax, length.out = 20), include.lowest = TRUE)
    data[,xfocus] <- sapply(data[,xfocus], get_midpoint)
    
    x_new <- data.frame(id = 1:length(unique(data[,xfocus])), variable = sort(unique(data[,xfocus])))
    Xbeta <- matrix(NA, nrow = n_iterations, ncol = nrow(x_new))
    for (i in 1:n_iterations) {
      Xbeta[i,] <- as.matrix(x_new$variable) %*% filter(coefs, .data$iteration == i)$value
    }
    coefs <- melt(Xbeta, varnames = c("iteration", "id")) %>%
      left_join(x_new, by = "id") %>%
      select(-id)    
  } else {
    data <- fit$data %>%
      mutate_at(vars(matches(xfocus)), factor)
  }
  
  # Influence
  influ <- get_influ(fit = fit, group = c(yfocus, xfocus), hurdle = hurdle)
  
  if (nrow(fit$ranef) > 0) {
    ylab1 <- "Coefficient"
  } else {
    ylab1 <- "Relative coefficient"
  }
  
  # Extract the legend on its own
  g2 <- function(a.gplot) {
    if (!is.gtable(a.gplot))
      a.gplot <- ggplotGrob(a.gplot)
    gtable_filter(a.gplot, 'guide-box', fixed = TRUE)
  }
  
  # Build the plot
  sp <- 0.05
  
  # The coefficients (top-left)
  p1 <- ggplot(data = coefs, aes(x = factor(.data$variable), y = exp(.data$value))) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_violin(colour = colour, fill = colour, alpha = 0.5, draw_quantiles = 0.5, scale = "width") +
    labs(x = NULL, y = ylab1) +
    scale_x_discrete(position = "top") +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = margin(b = sp, r = sp, unit = "cm"))
  
  # The bubble plot (bottom-left) and the legend for the bubble plot (top-right)
  p3a <- plot_bubble(df = data, group = c(yfocus, xfocus), sum_by = "row", xlab = xlab, ylab = ylab, zlab = "", fill = colour)
  p2 <- g2(p3a)
  p3 <- p3a +
    theme(legend.position = "none", plot.margin = margin(t = sp, r = sp, unit = "cm"), axis.text.x = element_text(angle = 45, hjust = 1))
  
  # The influence plot (bottom-right)
  p4 <- ggplot(data = influ, aes_string(x = as.character(yfocus))) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    # geom_violin(aes(y = .data$delta), colour = colour, fill = colour, alpha = 0.5, draw_quantiles = 0.5, scale = "width") +
    geom_violin(aes(y = exp(.data$delta)), colour = colour, fill = colour, alpha = 0.5, draw_quantiles = 0.5, scale = "width") +
    coord_flip() +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = "Influence") +
    theme_bw() +
    theme(legend.position = "none", plot.margin = margin(t = sp, l = sp, unit = "cm"))
  
  p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2, heights = c(1, 2), widths = c(2, 1))
}


################################################################################################
# Oxana;s draft
################################################################################################

#' Plot Coefficient Distribution and Influence (CDI)
#'
#' @description 
#' Generates a multi-pane layout visualising the coefficients, data distribution, 
#' and relative influence of a model predictor. 
#' 
#' @param fit A fitted model object (e.g., from `glm` or `gam`).
#' @param year Character string. The name of the temporal column (e.g., "fyear") 
#'   to be used as the focus variable on the y-axis of the distribution plot.
#' @param predictor 
#
#' 
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @importFrom stats as.formula coef
#' 
#' @export

plot_cdi <- function(fit, year = NULL, raw_data = NULL, predictor = NULL, compare_preds_list = NULL){
  
  # comparison switch
  
  compareOn <- !is.null(compare_preds_list)
  
  if (compareOn) compare_preds_df <- compare_preds_list$preds
  
  if (is.null(year)) {
    year <- get_first_term(fit = fit)
  }
  
  # Derive contribution of each term to response
  preds_list <- get_preds(fit, raw_data = raw_data)
  
  # Extract components from the preds list
  preds            <- preds_list$preds
  V                <- preds_list$V
  X_centered       <- preds_list$X_centered
  assigns          <- preds_list$assign
  all_model_terms  <- preds_list$terms
  
  # exclude year from terms
  terms_labels <- setdiff(all_model_terms, year)
  
  result <- setNames(lapply(terms_labels, function(term_label) {
    
    # Define col names for fitted terms and ses
    fit_colname <- sym(paste0('fit.', term_label))
    se_colname  <- sym(paste0('se.fit.', term_label))
    
    # Identify model matrix columns for this specific term
    match_idx <- which(all_model_terms == term_label)
    cols <- which(assigns == match_idx)
    
    # Extract the sub model matrix for this term
    Xi_centered <- X_centered[, cols, drop = FALSE]
    Vi <- V[cols, cols, drop = FALSE]
    
    # Define levels of term on which coefficient and distribution plots will be based
    # This is done by search for each column name in the data as a whole word in the
    # each term. This allows for matching of terms like 'poly(log(var),3)' with 'var'
    
    term_stripped <- all.vars(as.formula(paste("~", term_label))) 
    raw_values <- preds[[term_stripped]]
    
    # Numeric terms are cut into factors 
    if(is.numeric(raw_values)) {
      breaks <- pretty(raw_values, 30)
      step   <- diff(breaks[1:2])
      labels <- breaks + (step / 2)
      breaks <- c(breaks, max(breaks) + step)
      
      levels <- cut(raw_values, 
                    breaks = breaks, 
                    labels = labels, 
                    include.lowest = TRUE)
    } else {
      levels <- raw_values
    }
    
    coeffs <- preds %>%
      mutate(row_id = row_number()) %>%
      group_by(term = !!levels) %>%  
      summarise(
        coef = mean(!!fit_colname),
        se   = {
          # Calculate the average design matrix row for this bin
          X_bin_avg <- colMeans(Xi_centered[row_id, , drop = FALSE])
          # Project the covariance matrix: sqrt(avg_row %*% V %*% avg_row_T)
          sqrt(t(X_bin_avg) %*% Vi %*% X_bin_avg)
        },
        .groups = "drop"
      ) %>%
      mutate (
        lower = coef - (1.96 * se),
        upper = coef + (1.96 * se)
      )
    
    # Reorder levels according to coefficients for factors
    if(is.factor(raw_values)) {
      coeffs <- coeffs %>%
        arrange(coef) 
      
      coeffs$term <- factor(coeffs$term, levels = coeffs$term, ordered = TRUE)
      
      levels <- factor(levels, levels = coeffs$term, ordered = TRUE)
      
    }
    
    if(compareOn && term_stripped %in% names(compare_preds_df)){
      
      # Extract components from the comparison list
      comp_raw_values    <- compare_preds_df[[term_stripped]]
      comp_V             <- compare_preds_list$V
      comp_X_centered    <- compare_preds_list$X_centered
      comp_assign        <- compare_preds_list$assign
      comp_terms         <- compare_preds_list$terms
      
      # Numeric terms are cut into factors 
      if(is.numeric(comp_raw_values)) {
        
        comp_levels <- cut(comp_raw_values, 
                           breaks = breaks, 
                           labels = labels, 
                           include.lowest = TRUE)
      } else {
        comp_levels <- comp_raw_values
      }
      
      # Identify columns in the comparison model matrix
      comp_match_idx <- which(comp_terms == term_label)
      comp_cols <- which(comp_assign == comp_match_idx)
      # Extract the sub model matrix for this term
      comp_Xi_centered <- comp_X_centered[, comp_cols, drop = FALSE]
      comp_Vi <- comp_V[comp_cols, comp_cols, drop = FALSE]
      
      
      comp_coeffs <- compare_preds_df %>%
        mutate(row_id = row_number()) %>%
        group_by(term = !!comp_levels) %>%  
        summarise(
          coef = mean(!!fit_colname),
          se = {
            X_bin_avg <- colMeans(comp_Xi_centered[row_id, , drop = FALSE])
            as.numeric(sqrt(t(X_bin_avg) %*% comp_Vi %*% X_bin_avg))
          },
          .groups = "drop"
        ) %>%
        mutate (
          term = factor(term, levels = levels(coeffs$term), ordered = TRUE),
          lower = coef - (1.96 * se),
          upper = coef + (1.96 * se)
        )
    }
    
    # 1. Coefficients  Plot
    #_____________________________________________________________________________
    
    # Create plot theme depending on the term:
    
    if ((grepl('month|area|cluster', term_label, ignore.case=TRUE)) | (length(levels(levels)) > 8)) {
      
      # Vertical labels and tighter margins
      dynamic_theme <- theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
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
    x_labels <- if(grepl('vessel_key', term_label, ignore.case=TRUE)){
      NULL
      # }  else if(length(levels(levels)) > 15) {
      #   # Downsample labels for continuous vars otherwise they overlap, the logic is to keep every n-th label where n = original no of levels/10
      #   function(x) {
      #     x[seq_along(x) %% max(1, round(length(x) / 10)) != 1] <- ""
      #     return(x)
      #   }
    } else waiver()
    
    
    p1 <- ggplot(coeffs, aes(x = term, y = exp(coef))) +
      geom_point(shape = 2, size = 2) +
      geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), width = 0) +
      # Bottom Cap
      geom_segment(aes(x = as.numeric(term) - 0.05, xend = as.numeric(term) + 0.05, 
                       y = exp(lower), yend = exp(lower))) +
      # Top Cap
      geom_segment(aes(x = as.numeric(term) - 0.05, xend = as.numeric(term) + 0.05, 
                       y = exp(upper), yend = exp(upper))) +
      geom_hline(yintercept = 1, linetype = "dashed")+
      
      # Conditional block to display coeffs for a model to be compared with (e.g. last year)
      { if (compareOn && term_stripped %in% names(compare_preds_df)) {
        list(
          geom_point(data = comp_coeffs,
                   aes(x = term, y=exp(coef)),
                   shape = 17, size = 2, color = 'palevioletred4'),
          geom_errorbar(data = comp_coeffs, aes(ymin = exp(lower), ymax = exp(upper)), color = 'palevioletred4', width = 0) 
        )
      }} +
      
      scale_y_log10() +
      scale_x_discrete(drop = FALSE,
                       labels = x_labels, 
                       position = 'top') +
      theme_bw() +
      background_grid(major = "x", minor = "none") +
      labs(y = 'Coefficient', x = NULL) +
      dynamic_theme +
      theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
      )
    
    # 2. Distribution Plot
    #_____________________________________________________________________________
    
    distrs <- preds %>%
      group_by(focus = .[[year]], term = levels)  %>%
      summarise(count = n(), .groups = "drop_last") %>%
      
      mutate(
        total = sum(count),
        prop  = count / total
      ) %>%
      
      ungroup()
    
    p2 <- ggplot(distrs, aes(x = term, y = focus, size = sqrt(prop) * 20)) +
      geom_point(pch = 1) +
      scale_size_identity() +
      scale_x_discrete(labels = x_labels) +
      theme_bw() +
      background_grid(major = "xy", minor = "none") +
      labs(x = term_stripped, y = year) +
      dynamic_theme +
      theme(axis.line = element_line(colour = "black", linewidth = 0.8))
    
    
    # 3. Influence
    #______________________________________________________________________________
    
    # This is formula 3b from Bentley et al. 2012
    infl <- preds %>%
      group_by(level = !!sym(year)) %>%
      summarise(Influence = mean(!!sym(fit_colname))) 
    
    # Mean absolute magnitude of the term's impact, converted to annual % change
    overall <- exp(mean(abs(infl$Influence))) - 1
    
    # Annual growth % of this term's contribution
    n <- nrow(infl)
    trend <- exp(cov(1:n, infl$Influence) / var(1:n)) - 1
    
    infl %<>% 
      mutate(Influence = exp(Influence),
             Term = term_label, 
             overall = overall, 
             trend = trend) %>%
      relocate(Term, .after = level)
    
    if (compareOn && term_stripped %in% names(compare_preds_df)){
      
      comp_infl <- compare_preds_df %>%
        group_by(level = !!sym(year)) %>%
        summarise(Influence = exp(mean(!!sym(fit_colname)))) 
      comp_x_min <- min(comp_infl$Influence)
      comp_x_max <- max(comp_infl$Influence)
    }
    
    p3 <- ggplot(infl, aes(x = Influence, y = level)) +
      geom_hline(aes(yintercept = level), color = "grey", linewidth = 0.5) +
      geom_vline(xintercept = 1, linetype = "dashed") +
      
      # Conditional block to display coeffs for a model to be compared with (e.g. last year)
      { if (compareOn && term_stripped %in% names(compare_preds_df)) {
        list(
          # Outer shaded rectangle (+/- 0.05 padding)
          annotate("rect", 
                   xmin = comp_x_min - 0.05, xmax = comp_x_max + 0.05, 
                   ymin = -Inf, ymax = Inf, 
                   fill = "grey90", alpha = 0.5),
          
          # Inner shaded rectangle (exact range)
          annotate("rect", 
                   xmin = comp_x_min, xmax = comp_x_max, 
                   ymin = -Inf, ymax = Inf, 
                   fill = "grey80", alpha = 0.5),
          
          geom_point(data = comp_infl,
                     shape = 16, size = 3, color = 'palevioletred4'),
          geom_line(data = comp_infl, group = 1, linetype = 'dashed', color = 'palevioletred4')  
        )
      }} +
      
      geom_line(group = 1) +                    # Connects the dots in the order of the data
      geom_point(size = 3, pch = 16) +
      scale_x_continuous(labels = scales::label_number(accuracy = 0.01), expand = expansion(mult = c(0.2, 0.2))) +
      scale_y_discrete( position = "right") +
      labs(x = "Influence", y = NULL)+
      theme_bw() +
      background_grid(major = "y", minor = "none") +
      theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
      )
    
    
    pv <- ggplot() + theme_void()
    
    # Legend for comparison plot
    
    if(compareOn && term_stripped %in% names(compare_preds_df)) {
      
      labels = c("Current", 
                 paste("Previous update to", max(as.numeric(as.character(comp_infl$level)))),
                 "Range of previous update", 
                 "Range of \n previous update +5%")
      
      pv <- ggplot() +
        # Row 1: Current
        annotate("segment", x = 5, xend = 15, y = 50, yend = 50, color = "black", linetype = "solid") +
        annotate("point", x = 10, y = 50, shape = 16, size = 2, color = "black") +
        annotate("point", x = 3, y = 50, shape = 2, size = 2, color = "black") +
        
        # Row 2: Previous
        annotate("segment", x = 5, xend = 15, y = 40, yend = 40, color = "brown", linetype = "dotted") +
        annotate("point", x = 10, y = 40, shape = 16, size = 2, color = "brown") +
        annotate("point", x = 3, y = 40, shape = 17, size = 2, color = "brown") +
        
        # Row 3 & 4: Ranges
        annotate("point", x = 10, y = 30, shape = 15, size = 6, color = "grey80", alpha = 0.5) +
        annotate("point", x = 10, y = 20, shape = 15, size = 6, color = "grey90", alpha = 0.5) +
        
        # Text Labels
        annotate("text", x = 16, y = seq(from = 50, to = 20, by = -10), label = labels, hjust = 0, size = 4) +
        
        
        coord_cartesian(xlim = c(0, 60), ylim = c(0, 60), expand = FALSE) +
        
        theme_void() 
    }
    
    combined_plot <- p1 + pv + p2 + p3 + plot_layout(nrow = 2, ncol = 2, heights = c(1, 2), widths = c(2, 1)) &
      theme(
        # axis.title = element_text(size = 14), 
        #  axis.text = element_text(size = 14),  
        plot.margin = margin(0, 0, 0, 0)  
      )
    
    # print(combined_plot)
    
    
  }), terms_labels)
  
  return (result)
}
