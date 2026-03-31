#' Plot Residual Implied Index
#'
#' @param fit A fitted model object
#' @param grouping_var A character string specifying the column name to group by. 
#' Defaults to 'stat_area'.
#' @importFrom splines ns
#' @importFrom cowplot theme_cowplot
#' @return A ggplot2 object showing the implied vs. scaled indices.
#' @export


plot_RIC <- function(fit, grouping_var = 'stat_area', min_records = 10,  add.rho = TRUE){
  year <- get_first_term(fit)
  grouping_var <- sym(grouping_var)
  
  component <- if(length(fit$family$family)==2) {
    'Combined'
  } else if (!is.null(fit$family$family) && (fit$family$family %in% c("bernoulli", "binomial"))){
    'Binomial'
  } else {
    'Positive'
  }
  
  
  idx <- get_index(fit, format = "wide")
  idx <- idx %>%
    filter(Index==component)
  
  
  raw_data <- fit$data
  
  
  if (is.null(raw_data) && inherits(fit, "survreg")) {
    raw_data <- eval(fit$call$data)
  }
  
  ric_data <- tibble(level = raw_data[[year]],
                     !!grouping_var := raw_data[[as.character(grouping_var)]],                              
                     resid  = residuals(fit, type = 'response')) %>%
    add_count(level, !!grouping_var) %>% 
    filter(n >= min_records) %>%
    left_join(idx, by = 'level') %>%
    mutate(implied = if (grepl("log", formula(fit)[2])) {
      exp(resid + log(stan_unscaled))
    } else {
      resid + stan_unscaled
    }) %>%
    group_by(!!grouping_var) %>%    
    mutate(base_imp = {
      # Calculate arithmetic mean for each year 
      yearly_means <- tapply(implied, level, mean, na.rm = TRUE)
      # Calculate geometric mean of those yearly means
      exp(mean(log(yearly_means), na.rm = TRUE))
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
    mutate(level  = as.integer(as.character(level))) %>%
    complete(level = full_seq(level, 1)) %>% 
    ungroup()
    
  imp_count <- ric_data %>%
    group_by(!!grouping_var) %>%
    summarise(Num = sum(n, na.rm = T),
           rho=cor(implied, idx_scaled, use="pairwise.complete.obs"))
  
  
  p <-   ggplot(ric_data,
                aes(x=level,
                    y=imp_scaled))+
    geom_point(aes(size=n), alpha = 0.5)+
    geom_line()+
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
    scale_x_continuous(breaks = function(x) {
      # x is the range of years present in the data
      years <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
      
      if (length(years) > 15) {
        # If more than 15 years, show every 2nd year
        round(seq(min(years), max(years), by = 2))
      } else {
        # Otherwise show all
        round(years)
      }
    }) +
    theme_cowplot()+
    theme(axis.text.x = element_text(hjust = 0,
                                     angle = 90, size = 10)) +
    (if(add.rho) {
      list(
          geom_text(data = imp_count, aes(x = Inf, y = Inf, label = paste0("N = ", Num)), 
                    vjust = 1.2, hjust = 1.1, colour = 'deepskyblue4'),
          geom_text(data = imp_count, aes(x = Inf, y = Inf, label = paste0('rho == ', round(rho, 2))), 
                    vjust = 2.6, hjust = 1.1, colour = 'deepskyblue4', parse = TRUE)
        )
    } else {
      NULL
    })
  
  return(p)
}

