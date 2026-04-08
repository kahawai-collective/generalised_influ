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

plot_DHARMares <- function(fit){  
  #______________________________________________________________________________
  #  Prepare residuals 
  #______________________________________________________________________________
  
  diag_metrics <- data.frame(
    # the reason for sub-setting to [1:nrow(data)] is to accommodate survreg which has extra column in response
    observed = as.vector(model.response(model.frame(fit)))[1:nrow(model.frame(fit))],
    fitted   = predict(fit, type = "response")
  )
  n <- nrow(diag_metrics)
  if  (inherits(fit, "survreg")) {
    
    simulatedResponse <- replicate(1000, rweibull(nrow(diag_metrics), 
                                                  shape = 1/fit$scale, 
                                                  scale = predict(fit, type = "lp") %>% exp()))
    
    DHARMa <- DHARMa::createDHARMa(simulatedResponse = simulatedResponse, 
                                   observedResponse = diag_metrics$observed,
                                   fittedPredictedResponse = diag_metrics$fitted,
                                   seed = 123)
    diag_metrics$dharma_res <- DHARMa$scaledResiduals
    
  } else diag_metrics$dharma_res <- residuals(DHARMa::simulateResiduals(fit, n = 1000, refit = FALSE, plot = FALSE, seed = 123))
  
  # Transform fitted values to ranks
  diag_metrics <- diag_metrics %>%
    mutate(rank_fitted =  rank(fitted, ties.method = "average")/n,
           is_outlier = factor(dharma_res == 0 | dharma_res == 1, 
                               levels = c(FALSE, TRUE)))
  
  # prepare quantile test and colour for quantile regression line
  q_test <-  DHARMa::testQuantiles(diag_metrics$dharma_res, diag_metrics$rank_fitted, plot = F)
  p_vals <- q_test$pvals
  line_colours <- ifelse(p_vals <= 0.05, "red", "black")
  
  #______________________________________________________________________________
  #  Plotting code
  #______________________________________________________________________________
  
  # (1) ---Histogram of residuals versus uniform distribution---
  d1 <- ggplot(diag_metrics) + geom_histogram(aes(x=qnorm(dharma_res),
                                                  y = after_stat(density)),
                                              stat='bin',
                                              binwidth=0.05,
                                              fill =NA,
                                              colour='black') +
    labs(x='DHARMa residual',
         y='Density') +
    stat_function(fun = dnorm, 
                  #     args = list(min = 0, max = 1), 
                  color = "red", linewidth = 1)+
    theme_classic() +
    theme(axis.title=element_text(size=10))
  
  # (2) ---Residuals vs fitted values---
  
  # built-in function from DHARMa package, but result cannot be stored and used in patchwork
  # d2 <- DHARMa::plotResiduals(DHARMa, main = NULL)
  
  d2 <- ggplot(diag_metrics, aes(x = rank_fitted, y = dharma_res)) + 
    geom_point(aes(shape = is_outlier, 
                   color = is_outlier, 
                   alpha = is_outlier),
               fill = 'white',
               stat = 'unique') +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 8)) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    scale_alpha_manual(values = c("FALSE" = 0.2, "TRUE" = 1)) +
    labs(x='Model predictions (rank transformed)',
         y='DHARMa residual') +
    geom_hline(yintercept = c(0.25, 0.5, 0.75), color = line_colours, linetype = "dashed", alpha = 0.5) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
                color = line_colours[2], se = FALSE) +
    guides(shape = "none", color = "none", alpha = "none")+
    theme_classic() +
    theme(axis.title=element_text(size=10))
  
  # (3) ---QQ plot---
  
  # built-in function from DHARMa package, but result cannot be stored and used in patchwork
  # d3 <- DHARMa::plotQQunif(DHARMa, testDispersion = FALSE,
  #                          testUniformity = FALSE,
  #                          testOutliers = FALSE, main = NULL)
  
  d3 <- ggplot(diag_metrics) + geom_qq(aes(sample=dharma_res), distribution = stats::qunif, dparams = list(min = 0, max = 1),
                                       cex=0.3) +
    labs(x='DHARMa resid. theoretical quantile',
         y='DHARMa resid. sample quantile') +
    geom_abline(intercept=0,
                linetype='dotted',
                colour='blue') +
    theme_classic() +
    theme(axis.title=element_text(size=10))
  
  # (4) ---Observed versus fitted values--- 
  d4 <- ggplot(diag_metrics) + geom_point(aes(fitted,
                                              observed),
                                          cex=0.3,
                                          shape = 21,
                                          color = 'black',
                                          fill = 'white',
                                          stat = 'unique')  +
    labs(x='Fitted value',
         y='Observed value') +
    geom_abline(intercept=0,
                linetype='dotted',
                colour='blue') +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') +
    theme_classic() +
    theme(axis.title=element_text(size=10))
  
  p <-  plot_grid(d1, d2, d3, d4, nrow=2, align='hv')
  
  return(p)
}