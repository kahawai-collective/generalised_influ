#' Get the unstandardised indices
#' 
#' @param fit An object of class \code{brmsfit}.
#' @param year The year or time label (e.g. year, Year, fishing_year, etc).
#' @param rescale How to re-scale the series. Choose from "raw" to retain the raw unstandardised series, or a number to re-scale by. 
#' @param predictor =1 or =2 if a hurdle model is used. If NULL, combined index will be returned.
#' @return a \code{data.frame} or a \code{ggplot} object.
#' @importFrom brms is.brmsfit
#' @import dplyr
#' @export
#' 
get_unstandardised <- function(fit, year = NULL, rescale = 1, predictor = NULL) {
  
  if  (!any(class(fit) %in% c("sdmTMB", "glm", "brmsfit", "survreg"))) stop("This model class is not supported.")
  
  is_sdm <- inherits(fit, 'sdmTMB')
  
  if (is.null(year)) {
    year <- get_first_term(fit = fit)
  }
  
  if(inherits(fit, 'brmsfit')){
    
    response_name <- formula(fit)$formula[[2]]  # 2 is a position of response variable in formulas
    observed <- fit$data[[response_name]]
    
  } else if (is_sdm){
    
    observed <- fit$response[,1]
    observed_pos <- fit$response[,2]
    
  } else if (inherits(fit, "glm")){
    
    response_name <- as.character(formula(fit)[2])
    observed <- fit$model[[response_name]]
    
  } else if (inherits(fit, "survreg")){
    
    response_name <- as.character(formula(fit)[2])
    observed <- fit$model[[response_name]]
    observed <- as.numeric(observed[,1])
    
  }
  
  year_vec <- model.frame(fit)[[year]]
  
  
  
  # (1) Binomial component only
  if (!is.null(fit$family$family) && any(fit$family$family %in% c("bernoulli", "binomial") | grepl("hurdle", fit$family$family))) {
    
    indices_bin <- aggregate(list(unstan_prob=observed),                        # derive mean by year
                             list(level=year_vec),
                             mean) %>%
      mutate(unstan = exp(log(unstan_prob)-mean(log(unstan_prob))))         # derive ratio of each index to its geo mean (relative index)
    
    indices <-indices_bin %>% rename(unstan_unscaled = unstan_prob)
    
    # (2) sdmTMB model: Add positive component and combine  
    if(is_sdm){
      indices_pos <- aggregate(list(unstan=log(observed_pos)),                        # derive mean by year
                               list(level=year_vec),
                               mean, 
                               na.rm = TRUE) %>%                                    #removing NA values in positive model
        mutate(unstan = exp(unstan-mean(unstan)))
      indices_comb <- left_join(indices_bin, indices_pos, by = "level", suffix = c("_bin", "_pos")) %>%
        mutate(unstan_combined = unstan_bin * unstan_pos)           
      #assign binomial, positive or combined index depending on 'predictor' param
      indices <- if (is.null(predictor)) {
        indices_comb
      } else {
        switch(as.character(predictor),
               "1" = indices_bin,
               "2" = indices_pos,
               stop("Invalid predictor value")
        )
      }
    }
    
    
    # Darcy's code
    # prop <- data.frame(y = fit$data[,1], Year = fit$data[,year]) %>%
    #   mutate(y = ifelse(y > 0, 1, 0)) %>%
    #   group_by(Year) %>%
    #   summarise(p = sum(y) / n())
    
    # unstd <- data.frame(y = fit$data[,1], Year = fit$data[,year]) %>%
    #   filter(y > 0) %>%
    #   group_by(Year) %>%
    #   summarise(cpue = exp(mean(log(y)))) %>%
    #   left_join(prop, by = "Year") %>%
    #   mutate(cpue = cpue * p)
    
  } else {
    
    # (3) Positive component only
    logged <- any(grepl('log', as.character(response_name), fixed = TRUE))

    if(logged) log_observed = observed else log_observed = log(observed)
    
    indices <- aggregate(list(unstan_pos=log_observed),                        # derive mean by year
                         list(level=year_vec),
                         mean) %>%
      mutate(unstan = exp(unstan_pos-mean(unstan_pos)))  %>%             # derive ratio of each index to its geo mean (relative index)
      select(-unstan_pos)
    # Darcy's code
    # unstd <- data.frame(y = fit$data[,1], Year = fit$data[,year]) %>%
    #   group_by(Year) %>%
    #   summarise(cpue = exp(mean(log(y))))
  }
  
  
  return(indices)
}
# Darcy's code
#   gm <- geo_mean(unstd$cpue)
# 
#   fout <- unstd %>%
#     mutate(Mean = cpue, Median = cpue) %>%
#     select(-cpue)
# 
#   # Rescale the series
#   if (rescale == "raw") {
#     # nothing to do
#   } else if (is.numeric(rescale)) {
#     fout$Mean <- fout$Mean / gm * rescale
#     fout$Median <- fout$Median / gm * rescale
#   }
# 
#   return(fout)
# }


#' Get the standardised indices
#' 
#' Get the standardised indices each year with associated uncertainty and return either as a table or a ggplot.
#' 
#' @param fit An object of class \code{brmsfit}.
#' @param year The year or time label (e.g. year, Year, fishing_year, etc).
#' @param probs The percentiles to be computed by the \code{quantile} function.
#' @param rescale How to re-scale the series. Choose from "raw" to retain the raw series, "unstandardised" to re-scale to the geometric mean of the unstandardised series, or a number to re-scale by. 
#' @param do_plot Return a \code{ggplot} object instead of a \code{data.frame}.
#' @param ... Additional parameters passed to \code{fitted}.
#' @return a \code{data.frame} or a \code{ggplot} object.
#' @importFrom stats fitted
#' @importFrom brms is.brmsfit
#' @import insight
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @export
#' 
get_index <- function(fit, year = NULL, probs = c(0.025, 0.975), rescale = 1, pred_grid = NULL, predictor = NULL, ...) {
  
  if  (!inherits(fit, c("sdmTMB", "glm", "survreg", "brmsfit"))) stop("This model class is not supported.")
  
  ##################################
  # Get all the necessary variables
  ##################################
  is_sdm <- inherits(fit, 'sdmTMB')
  
  # Unstan indices
  indices <- get_unstandardised(fit = fit, predictor = predictor)
  
  # Name if first term
  if (is.null(year)) {
    year <- get_first_term(fit = fit)
  }
  # levels of fyear
  yrs <-  sort(unique(model.frame(fit)[[year]]))
  
  # number of years
  n <- length(yrs)
  
  # model data
  raw_data <- insight::get_data(fit)
  
  # name of the response variable
  response_name <- insight::find_response(fit)
  
  ###################################
  # Create  data for prediction
  ###################################
  
  newdata <- raw_data %>%
    dplyr::select(-all_of(year), -all_of(response_name)) %>%
    summarise(across(everything(), ~mean_or_mode(.x))) %>%
    tidyr::expand_grid(!!year := yrs )
  
  
  if (is.brmsfit(fit)) {
    
    draws <- fitted(object = fit, newdata = newdata, probs = c(probs[1], 0.5, probs[2]), re_formula = NA, scale = "response", summary = F)
    
    colnames(draws) <- yrs
    
    draws %<>%
      as.data.frame() %>%
      mutate(.iteration = row_number()) %>%
      pivot_longer(cols = -.iteration,
                   names_to = 'level',
                   values_to = '.value') 
    
    # Get the predicted CPUE by year
    # fout1 <- fitted(object = fit, newdata = newdata, probs = c(probs[1], 0.5, probs[2]), re_formula = NA) %>% 
    #   data.frame() %>%
    #   rename(Qlower = 3, Qupper = 5) %>% # this renames the 3rd and the 5th columns
    #   mutate(CV = Est.Error / Estimate) %>% # CV = SD / mu
    #   mutate(Year = yrs) %>%
    #   mutate(Model = as.character(fit$formula)[1], Distribution = as.character(fit$family)[1], Link = as.character(fit$family)[2])
    
    # # Rescale the predicted CPUE. The options are:
    # # 1. raw - don't rescale
    # # 2. one - rescale so that the series has a geometric mean of one
    # # 3. unstandardised - rescale to the geometric mean of the unstandardised series
    # # 4. a user defined number
    # fout <- fout1
    # if (rescale == "raw") {
    #   # nothing to do
    # } else if (rescale == "unstandardised") {
    #   unstd <- get_unstandarsied(fit = fit, year = year, rescale = "raw")
    #   gm <- geo_mean(unstd$Mean)
    #   fout$Estimate <- fout$Estimate / geo_mean(fout$Estimate) * gm
    #   fout$Qlower <- fout$Qlower / geo_mean(fout$Q50) * gm
    #   fout$Qupper <- fout$Qupper / geo_mean(fout$Q50) * gm
    #   fout$Q50 <- fout$Q50 / geo_mean(fout$Q50) * gm
    # } else if (is.numeric(rescale)) {
    #   fout$Estimate <- fout$Estimate / geo_mean(fout$Estimate) * rescale
    #   fout$Qlower <- fout$Qlower / geo_mean(fout$Q50) * rescale
    #   fout$Qupper <- fout$Qupper / geo_mean(fout$Q50) * rescale
    #   fout$Q50 <- fout$Q50 / geo_mean(fout$Q50) * rescale
    # }
    # fout$Est.Error <- fout$CV * fout$Estimate # SD = CV * mu
    
    # if (do_plot) {
    #   p1 <- ggplot(data = fout1, aes(x = Year)) +
    #     geom_errorbar(aes(y = Q50, ymin = Qlower, ymax = Qupper)) +
    #     geom_point(aes(y = Q50)) +
    #     geom_errorbar(aes(y = Estimate, ymin = Estimate - Est.Error, ymax = Estimate + Est.Error), colour = "red", alpha = 0.75) +
    #     geom_point(aes(y = Estimate), colour = "red", alpha = 0.75) +
    #     theme_bw()
    #   
    #   p2 <- ggplot(fout, aes(x = Year)) +
    #     geom_errorbar(aes(y = Q50, ymin = Qlower, ymax = Qupper)) +
    #     geom_point(aes(y = Q50)) +
    #     geom_errorbar(aes(y = Estimate, ymin = Estimate - Est.Error, ymax = Estimate + Est.Error), colour = "red", alpha = 0.75) +
    #     geom_point(aes(y = Estimate), colour = "red", alpha = 0.75) +
    #     theme_bw()
    #   return(p1 + p2)
    # } else {
    # Rename and reorder columns
    # index_brms <- fout %>%
    #   rename(Mean = Estimate, SD = Est.Error, Median = Q50) %>%
    #   mutate(Qlow = Qlower, Qup = Qupper) %>%
    #   rename_with(~paste0("Q", probs[1] * 100), Qlow) %>%
    #   rename_with(~paste0("Q", probs[2] * 100), Qup) %>%
    #   relocate(Year, Mean, SD, CV)
    #   } 
    
    
    
  }else if(inherits(fit, c("glm", "survreg"))){
    # GLMs and Survreg. 
    
    # Extract model coeffs
    cfs <- coefficients(fit)
    
    # Extract covariance matrix. Survreg has extra param Log(scale), we subset vcov matrix to the size fo coefficients only 
    V <- vcov(fit)[1:length(cfs), 1:length(cfs)]
    
    # Create draws of model coeffs
    beta_draws <- mvtnorm::rmvnorm(1000,cfs,V)
    
    # Force the levels in newdata to match the levels in raw_data
    for (col_name in names(newdata)) {
      if (is.factor(raw_data[[col_name]])) {
        newdata[[col_name]] <- factor(newdata[[col_name]], 
                                      levels = levels(raw_data[[col_name]]))
      }
    }
    
    # Model matrix for new data
    X <- model.matrix(delete.response(terms(fit)), data = newdata)
    
    # Predicted response for new data
    draws <- exp(beta_draws %*% t(X))
    
    # Different process for binomial component
    if (!is.null(fit$family$family) && any(fit$family$family %in% c("bernoulli", "binomial"))){
      
      rows = c(1,which(substr(row.names(V),
                              1,
                              nchar(year))==year))
      
      resp_name <- names(model.frame(fit))[1]
      int <- logit(aggregate(list(unstan=fit$model[[resp_name]]),
                             list(level=fit$model[[year]]),
                             mean)$unstan)[1]
      cfs[1] <- int
      
      bin <- mvtnorm::rmvnorm(1000,cfs,V)[,rows]
      bin <- cbind(inv_logit(bin[,1]),inv_logit(bin[,1] + bin[,2:length(yrs)]))
      
      draws <- bin
    }
    
    # Rearrange df for index derivation
    colnames(draws) <- yrs
    
    draws %<>%
      as.data.frame() %>%
      mutate(.iteration = row_number()) %>%
      pivot_longer(cols = -.iteration,
                   names_to = 'level',
                   values_to = '.value') 
    
  }else if(is_sdm) {
    # Spatio-temporal models 
    
    # if (is.null(predictor)) {
    #   # Combined index
    #   
    #   predict_both <- predict(spatiotemporal, newdata = pred_grid, return_tmb_object = TRUE)
    #   index_sdm <- sdmTMB::get_index(predict_both,  bias_correct = TRUE) %>%
    #     rename(level = !!sym(year)) %>%
    #     mutate(stan = exp(log_est - mean (log_est)),
    #            stanLower = exp(log_est - mean (log_est) - 1.96*se),
    #            stanUpper = exp(log_est - mean (log_est) + 1.96*se)
    #     )%>%
    #     select(level, stan, stanLower, stanUpper)
    #   
    #   
    # } else {
    #   # Non combined index
    
    if (predictor == 1){
      # Binomial Index
      
      predict_sim <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE, nsim = 1000, model = 1, type = "response")
      
    } else {
      # Combined index or Positive Index
      
      model <- if (length(predictor) > 0) 2 else NA
      predict_sim <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE, nsim = 1000, model = model)
      
    }
    
    # generate draws of indices
    draws <- sdmTMB::get_index_sims(predict_sim, return_sims = T) %>% # need to think about the area here
      rename(level = !!sym(year)) 
    
  }
  
  # This step is common to all models
  # Do geometric mean transformation to draws, then summarise to find medians and quantiles
  
  index_stan <- draws %>%
    
    group_by(.iteration) %>%
    mutate(
      level = factor(level),
      rel_idx = exp(log(.value) - mean(log(.value)))
    ) %>%
    ungroup() %>%
    group_by(level) %>%
    summarise(
      stan_unscaled = median(.value),
      stanLower_unscaled = quantile(.value, 0.025),
      stanUpper_unscaled = quantile(.value, 0.975),
      stan = median(rel_idx),
      stanLower = quantile(rel_idx, 0.025),
      stanUpper = quantile(rel_idx, 0.975)
    )
  
  indices <- indices %>%
    left_join(index_stan,  by = 'level')
  
  return(indices)
  
}

##################################################################################
# Draft approach for binomial index where index is predicted from the model at modes for factors 18.02.2026
##################################################################################

# # GLMs
# 
# # Extract covariance matrix
# V <- summary(fit)$cov.scaled
# 
# # Extract model coeffs
# cfs <- coefficients(fit)
# 
# # Create draws of model coeffs
# beta_draws <- mvtnorm::rmvnorm(1000,cfs,V)
# 
# # Force the levels in newdata to match the levels in raw_data
# for (col_name in names(newdata)) {
#   if (is.factor(raw_data[[col_name]])) {
#     newdata[[col_name]] <- factor(newdata[[col_name]], 
#                                   levels = levels(raw_data[[col_name]]))
#   }
# }
# 
# # Model matrix for new data
# X <- model.matrix(delete.response(terms(fit)), data = newdata)
# 
# # Predicted response for new data
# draws <- beta_draws %*% t(X)
# 
# # Bring to log-scale for binomial component
# if (!is.null(fit$family$family) && any(fit$family$family %in% c("bernoulli", "binomial"))){
#   draws <- log(inv_logit(draws))
# }
# 
# # Prepare df for index derivation
# colnames(draws) <- yrs
# 
# draws %<>%
#   as.data.frame() %>%
#   mutate(.iteration = row_number()) %>%
#   pivot_longer(cols = -.iteration,
#                names_to = 'level',
#                values_to = '.value') 
