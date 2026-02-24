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
#' @importFrom stats fitted
#' @import brms
#' @import ggplot2
#' @import dplyr
#' @export
#' 
plot_index <- function(index, 
                       year = NULL, 
                       fill = "black", 
                       probs = c(0.25, 0.75),
                       rescale = 1,
                       predictor = NULL,
                       show_unstandardised = TRUE) {
  
  if (is.null(year)) {
    year <- get_first_term(fit = fit)
  }
  
 index_long <- index %>%
    rename(Standardised = stan, Unstandardised = unstan) %>%
    pivot_longer(
      cols = c(Standardised, Unstandardised), 
      names_to = "index_type", 
      values_to = "Index"
    )
  
  if (!show_unstandardised) {
    index_long <- index_long %>% filter(index_type != "Unstandardised")
  }
  
 p <-  ggplot(index_long, aes(x = level, y = Index, group = index_type)) +
    geom_ribbon(aes(ymin = stanLower, ymax = stanUpper), 
                data = filter(index_long, index_type == "Standardised"),
                fill = fill, color = NA, alpha = 0.1) +
    geom_line(aes(linetype = index_type, colour = index_type)) +
    geom_point(aes(shape = index_type, colour = index_type), size = 3) +
    scale_shape_manual(values = c("Standardised" = 16, "Unstandardised" = 1))+
    scale_colour_manual(values = c( fill, "grey40")) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Fishing year", y = "Index") +
    theme_cowplot() +
    theme(
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.position = c(0.05, 0.95), 
      legend.justification = c(0, 1), 
      legend.title = element_blank(),
      legend.key.width = unit(2, "cm"),
      legend.background = element_rect(fill = "white", color = NA), 
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      
    )
  
  return(p)
}


#' Plot the hurdle and positive components
#' 
#' @param fit An object of class \code{brmsfit}.
#' @param year the year or time label.
#' @param fill the fill colour for the percentiles.
#' @param probs The percentiles to be computed by the \code{quantile} function.
#' @return a \code{ggplot} object.
#' 
#' @author Darcy Webber \email{darcy@quantifish.co.nz}
#' 
#' @importFrom stats fitted
#' @import brms
#' @import ggplot2
#' @import dplyr
#' @export
#' 
plot_hurdle <- function(fit, year = "Year", fill = "purple", probs = c(0.025, 0.975)) {
  
  if (!is.brmsfit(fit)) stop("fit is not an object of class brmsfit.")
  
  yrs <- sort(unique(fit$data[,year]))
  n <- length(yrs)
  
  # Create newdata for prediction (using fitted)
  newdata <- fit$data %>% slice(rep(1, n))
  for (j in 1:ncol(newdata)) {
    x <- fit$data[,j]
    newdata[,j] <- ifelse(is.numeric(x), mean(x), NA)
  }
  newdata[,year] <- yrs
  
  # Get the positive component
  mu <- fitted(object = fit, newdata = newdata, probs = c(probs[1], 0.5, probs[2]), re_formula = NA, dpar = "mu") %>% 
    cbind(newdata) %>%
    rename(Qlower = 3, Qupper = 5) %>% # this renames the 3rd and the 5th columns
    mutate(model = "Positive") %>%
    mutate(Estimate = exp(.data$Estimate), Q50 = exp(.data$Q50), Qlower = exp(.data$Qlower), Qupper = exp(.data$Qupper))
  # mu$Estimate <- mu$Estimate / geo_mean(mu$Estimate)
  mu$Qlower <- mu$Qlower / geo_mean(mu$Q50)
  mu$Qupper <- mu$Qupper / geo_mean(mu$Q50)
  mu$Q50 <- mu$Q50 / geo_mean(mu$Q50)

  # Get the hurdle component
  hu <- fitted(object = fit, newdata = newdata, probs = c(probs[1], 0.5, probs[2]), re_formula = NA, dpar = "hu") %>% 
    cbind(newdata) %>%
    rename(Qlower = 3, Qupper = 5) %>% # this renames the 3rd and the 5th columns
    mutate(model = "Hurdle")# %>%
    # mutate(Estimate = inv_logit(Estimate), Qlower = inv_logit(Qlower), Qupper = inv_logit(Qupper))
  # hu$Estimate <- hu$Estimate / geo_mean(hu$Estimate)
  hu$Qlower <- hu$Qlower / geo_mean(hu$Q50)
  hu$Qupper <- hu$Qupper / geo_mean(hu$Q50)
  hu$Q50 <- hu$Q50 / geo_mean(hu$Q50)
  
  # Get the combined series
  # bt <- fitted(object = fit, newdata = newdata, probs = c(probs[1], 0.5, probs[2]), re_formula = NA) %>%
  #   cbind(newdata) %>%
  #   rename(Qlower = 3, Qupper = 5) %>% # this renames the 3rd and the 5th columns
  #   mutate(model = "Combined")
  # bt$Qlower <- bt$Qlower / geo_mean(bt$Q50)
  # bt$Qupper <- bt$Qupper / geo_mean(bt$Q50)
  # bt$Q50 <- bt$Q50 / geo_mean(bt$Q50)
  bt <- get_index(fit = fit, year = year, probs = probs, rescale = 1) %>%
    mutate(model = "Combined")
  bt[year] <- bt$Year
  
  df <- bind_rows(mu, hu, bt)
  # df$model <- factor(df$model, levels = c("Hurdle", "Positive", "Combined"))
  
  p <- ggplot(data = df, aes(x = .data$year, y = .data$Q50, group = .data$model)) +
    geom_ribbon(aes(ymin = .data$Qlower, ymax = .data$Qupper, fill = .data$model), alpha = 0.5, colour = NA) +
    geom_line(aes(colour = .data$model, linetype = .data$model)) +
    geom_point(aes(colour = .data$model)) +
    labs(x = NULL, y = "Index") +
    # scale_colour_manual(values = c("grey", fill)) +
    # scale_fill_manual(values = c("grey", fill)) +
    # scale_linetype_manual(values = c("dashed", "solid")) +
    # scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    theme_bw() +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank(), legend.key.width = unit(2, "cm")) +
    guides(color = guide_legend(override.aes = list(fill = NA)))

    return(p)
}