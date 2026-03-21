#' Get the unstandardised indices
#' 
#' @param fit An object of class \code{brmsfit}.
#' @param year The year or time label (e.g. year, Year, fishing_year, etc).
#' @param rescale How to re-scale the series. Choose from "raw" to retain the raw unstandardised series, or a number to re-scale by. 
#' @param predictor =1 or =2 if a hurdle model is used. If NULL, combined index will be returned.
#' @return a \code{data.frame} or a \code{ggplot} object.
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
    
    
    } else {
    
    # (3) Positive component only
    logged <- any(grepl('log', as.character(response_name), fixed = TRUE))

    if(logged) log_observed = observed else log_observed = log(observed)
    
    indices <- aggregate(list(unstan_pos=log_observed),                        # derive mean by year
                         list(level=year_vec),
                         mean) %>%
      mutate(unstan = exp(unstan_pos-mean(unstan_pos)))  %>%             # derive ratio of each index to its geo mean (relative index)
      select(-unstan_pos)
    
  }
  
  
  return(indices)
}


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
  
  #_____________________________________________________________________________
  # Get all the necessary variables
  #_____________________________________________________________________________
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
  raw_data <- fit$data
  
  
  if (is.null(raw_data) && inherits(fit, "survreg")) {
    raw_data <- eval(fit$call$data)
  }
  
  # Forumula
  Formula <- if (is_sdm){
    
    if  (is.null(predictor)) stop("Argument 'predictor' is missing. Please specify 1 for the first part or 2 for the hurdle part.")
    formula(fit)[[predictor]]
    
  } else    { formula(fit)
    
  }
  
  mod_terms <- all.vars(as.formula(Formula))
  
  # name of the response variable
  response_name <- insight::find_response(fit)
  
  #_____________________________________________________________________________
  # Create  data for prediction
  #_____________________________________________________________________________
  
  
  # This bit not working, return to it later. Use dplyr for now.
  
  # cols_to_keep <- setdiff(mod_terms, c(year, response_name))
  # mod_data <- raw_data [, cols_to_keep, with = FALSE]
  # 
  # # Apply mean_or_mode to each column 
  # mean_mode_row <- as.data.frame(lapply(mod_data, mean_or_mode))
  # 
  # # Expand Grid: join with the 'yrs' vector
  # grid_list <- c(list(yrs), mean_mode_row)
  # names(grid_list)[1] <- year
  # newdata <- do.call(expand.grid, c(grid_list, stringsAsFactors = FALSE))
  # 
  
  cols_to_keep <- setdiff(mod_terms, c(year, response_name))
  
  # Subset the data to only terms that are in the model
  newdata <- raw_data %>%
    dplyr::select(all_of(cols_to_keep)) %>%
    summarise(across(everything(), ~mean_or_mode(.x))) %>%
    tidyr::expand_grid(!!year := yrs)
  
  #_____________________________________________________________________________
  # Draw from model predictions
  #_____________________________________________________________________________
  if (is.brmsfit(fit)) {
    
    draws <- fitted(object = fit, newdata = newdata, probs = c(probs[1], 0.5, probs[2]), re_formula = NA, scale = "response", summary = F)
    
    colnames(draws) <- yrs
    
    draws %<>%
      as.data.frame() %>%
      mutate(.iteration = row_number()) %>%
      pivot_longer(cols = -.iteration,
                   names_to = 'level',
                   values_to = '.value') 
    
    
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
        
        # Accommodate situation where factors were stores as characters, otherwise model matrix creation fails
      } else if (is.character(raw_data[[col_name]])){
        newdata[[col_name]] <- factor(newdata[[col_name]], 
                                      levels = sort(unique(as.character(raw_data[[col_name]]))))
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

#_____________________________________________________________________________
# Draft approach for binomial index where index is predicted from the model at modes for factors 18.02.2026
#_____________________________________________________________________________

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
