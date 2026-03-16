#' Plots predicted values versus residuals
#' 
#' This plots predicted values against residuals including uncertainty, using the fitted and residual functions on a brmsfit object.
#' 
#' @param fit An object of class \code{brmsfit}.
#' @param data a data frame with the same dimensions as the model data
#' @param year the year column
#' @param groups the grouping/facet variable
#' @return a ggplot object
#' 
#' @author Darcy Webber \email{darcy@quantifish.co.nz}
#' 
#' @importFrom stats fitted residuals
#' @import ggplot2
#' @import dplyr
#' @export
#' 
plot_implied_residuals <- function(fit, data = NULL, year = "Year", groups = "Species") {
  
  grp <- c(year, groups)
  grp2 <- c("Year", groups)
  
  # Get the data
  if (is.null(data)) {
    data <- fit$data
  }
  data <- data %>%
    mutate(id = 1:n()) %>%
    rename(Year = all_of(year))

  # Extract predicted values and normalise
  idx <- get_index(fit, year = year)
  idx$Median <- idx$Median - mean(idx$Median)

  # Extract residuals
  ires <- residuals(fit, summary = FALSE) %>% 
    melt() %>%
    rename(iteration = .data$Var1, id = .data$Var2, residual = .data$value) %>%
    left_join(data, by = "id") %>% 
    left_join(idx, by = "Year") %>%
    mutate(implied = .data$Median + .data$residual) %>%
    group_by_at(grp2) %>%
    summarise(Mean = mean(.data$implied), 
              Median = median(.data$implied),
              Qlower = quantile(.data$implied, probs = 0.05, na.rm = TRUE), 
              Qupper = quantile(.data$implied, probs = 0.95, na.rm = TRUE))

  p <- ggplot(data = ires, aes(x = .data$Year, y = .data$Median)) +
    geom_line(data = idx, aes(x = .data$Year, y = .data$Median), group = 1, linetype = "dashed") +
    geom_pointrange(aes(min = .data$Qlower, max = .data$Qupper), colour = "purple") +
    geom_line(group = 1, colour = "purple") +
    labs(x = NULL, y = "Residuals") +
    facet_wrap(as.formula(paste("~", groups)), ncol = 2) +
    theme_bw()
  
  return(p)
}



# plot_implied_residuals2 <- function(fit, data = NULL, year = "Year", groups = "Species") {
#   
#   grp <- c(year, groups)
#   grp2 <- c("Year", groups)
#   
#   # Get the data
#   if (is.null(data)) {
#     data <- fit$data
#   }
#   data <- data %>%
#     mutate(id = 1:n()) %>%
#     rename(Year = all_of(year))
#   
#   # Extract predicted values and normalise
#   idx <- get_index(fit, year = year)
#   idx$Estimate <- idx$Estimate - mean(idx$Estimate)
#   
#   # Extract residuals
#   ires <- residuals(fit, summary = FALSE) %>% 
#     melt() %>%
#     rename(iteration = .data$Var1, id = .data$Var2, residual = .data$value) %>%
#     left_join(data, by = "id") %>% 
#     left_join(idx, by = "Year") %>%
#     mutate(implied = .data$residual) %>%
#     group_by_at(grp2) %>%
#     summarise(Estimate = mean(.data$implied), 
#               Qlower = quantile(.data$implied, probs = 0.05, na.rm = TRUE), 
#               Qupper = quantile(.data$implied, probs = 0.95, na.rm = TRUE))
#   
#   p <- ggplot(data = ires, aes(x = .data$Year, y = .data$Estimate)) +
#     geom_pointrange(aes(min = .data$Qlower, max = .data$Qupper), colour = "purple") +
#     labs(x = NULL, y = "Residuals") +
#     facet_wrap(as.formula(paste("~", groups)), ncol = 2) +
#     theme_bw()
#   
#   return(p)
# }



#' Plots predicted values versus residuals
#' 
#' This plots predicted values against residuals including uncertainty, using the fitted and residual functions on a brmsfit object.
#' 
#' @param fit An object of class \code{brmsfit}.
#' @param trend show a loess smoother or linear line
#' @param type The type of the residuals, either "ordinary" or "pearson".
#' @return a ggplot object
#' 
#' @author Darcy Webber \email{darcy@quantifish.co.nz}
#' 
#' @importFrom stats fitted residuals
#' @import ggplot2
#' @import dplyr
#' @export
#' 
plot_predicted_residuals <- function(fit, trend = "loess", type = "ordinary") {
  # Extract predicted values
  pred <- fitted(fit) %>% 
    data.frame()
  names(pred) <- paste0("pred.", names(pred))
  
  # Extract residuals
  resid <- residuals(fit, type = type) %>% 
    data.frame()
  names(resid) <- paste0("resid.", names(resid))
  df <- cbind(resid, pred, fit$data)
  
  p <- ggplot(data = df, aes(x = .data$pred.Estimate, y = .data$resid.Estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmax = .data$pred.Q2.5, xmin = .data$pred.Q97.5, height = 0), alpha = 0.75) +
    geom_pointrange(aes(ymin = .data$resid.Q2.5, ymax = .data$resid.Q97.5), alpha = 0.75) +
    labs(x = "Predicted values", y = "Residuals") +
    theme_bw()
  
  if (trend == "loess") {
    p <- p + geom_smooth(method = "loess", se = FALSE, formula = y ~ x)
  }
  if (trend %in% c("lm", "linear")) {
    p <- p + geom_smooth(method = "lm", se = FALSE, formula = y ~ x)
  }  
  
  return(p)
}


#' Plot Residual Implied Index
#'
#' @param fit A fitted model object
#' @param grouping_var A character string specifying the column name to group by. 
#' Defaults to 'stat_area'.
#' @importFrom splines ns
#' @return A ggplot2 object showing the implied vs. scaled indices.
#' @export


plot_RIC <- function(fit, grouping_var = 'stat_area', add.rho = TRUE){
  year <- get_first_term(fit)
  grouping_var <- sym(grouping_var)
  idx <- get_index(fit)
  
  raw_data <- fit$data
  
  
  if (is.null(raw_data) && inherits(fit, "survreg")) {
    raw_data <- eval(fit$call$data)
  }
  
  ric_data <- tibble(level = raw_data[[year]],
                     !!grouping_var := raw_data[[as.character(grouping_var)]],                              
                      resid  = residuals(fit, type = 'response')) %>%
    left_join(idx, by = 'level') %>%
    mutate(implied = if (grepl("log", formula(fit)[2])) {
      exp(resid + log(stan_unscaled))
    } else {
      resid + stan_unscaled
    }) %>%
    group_by(!!grouping_var) %>%    
    mutate(base_imp = {
      # Calculate arithmetic mean for each year 
      yearly_means <- tapply(implied, level, mean)
      # Calculate geometric mean of those yearly means
      exp(mean(log(yearly_means)))
    },
           imp_scaled = implied / base_imp,
           idx_scaled = stan) %>%
     group_by(!!grouping_var, level) %>%
    summarise(n = n(),
              se = sd(imp_scaled)/ sqrt(n()),
              implied = mean(implied),
              idx = mean(stan_unscaled),
              imp_scaled = mean(imp_scaled),
              idx_scaled = mean(idx_scaled)) %>%
    mutate(Num = sum(n),
           rho=cor(implied, idx_scaled, use="pairwise.complete.obs"))
  
  p <-   ggplot(ric_data,
                aes(x=level,
                    y=imp_scaled))+
    geom_point(aes(size=n), alpha = 0.5)+
    geom_line(aes(group = 1))+
    geom_errorbar(aes(ymin=(imp_scaled-1.96*se),
                      ymax=(imp_scaled+1.96*se)),
                  size=0.3,
                  width=0.3)+
    geom_hline(yintercept=1,
               linetype=3,
               colour='grey')+
    geom_line(aes(y=idx_scaled, group = 1),
              col='grey')+
    facet_wrap(as.formula(paste("~", as.character(grouping_var))),
               ncol=2,scales='free_y')+
    labs(x='Fishing year',
         y='Index',
         size="Records") +
    theme_cowplot()+
    theme(axis.text.x = element_text(hjust = 0,
                                     angle = 90)) +
    (if(add.rho) {
      list(
        geom_text(aes(x = Inf, y = Inf, label = Num), 
                  vjust = 1.2, hjust = 1.1, colour = 'deepskyblue4'),
        geom_text(aes(x = Inf, y = Inf, label = paste0('rho == ', round(rho, 2))), 
                  vjust = 2.6, hjust = 1.1, colour = 'deepskyblue4', parse = TRUE)
      )
    } else {
      NULL
    })
  
  return(p)
}

