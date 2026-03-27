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
#' @importFrom splines ns
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



#________________________________________________________________________________________________________________
#' Extract and Combine Standardized Indices from Multiple Model Classes
#'
#' @description
#' A unified interface to extract temporal indices from \code{sdmTMB}, \code{glm}, 
#' \code{brmsfit}, and \code{survreg} models. This function specifically handles 
#' "Hurdle" or "Delta" structures by combining binomial (occurrence) and 
#' positive (magnitude) components into a single "Combined" expected value index, 
#' standardised by the geometric mean.
#'
#' @param fit A model object of class \code{sdmTMB}, \code{glm}, \code{brmsfit}, or \code{survreg}.
#' @param hurdle_fit Optional. A separate model object (e.g., a binomial GLM) if 
#'   the primary \code{fit} only represents the positive component of a delta model.
#' @param year Character string naming the temporal variable (e.g., "year"). 
#'   If \code{NULL}, the function attempts to automatically detect the first term.
#' @param probs Numeric vector of length 2 specifying the lower and upper quantiles 
#'   for uncertainty intervals. Defaults to \code{c(0.025, 0.975)}.
#' @param rescale Scale factor for the series. Defaults to 1.
#' @param pred_grid A \code{data.frame} for predictions. Required for \code{sdmTMB} 
#'   to define the spatial area for integration.
#' #'
#' @return A Dataframe containing indices
#' #' @importFrom stats coefficients vcov formula model.frame model.matrix delete.response terms aggregate quantile median
#' @importFrom dplyr select mutate group_by ungroup summarise across everything rename left_join bind_rows
#' @importFrom tidyr pivot_longer pivot_wider expand_grid
#' @importFrom rlang !! sym
#' @importFrom mvtnorm rmvnorm
#' @importFrom insight find_response
#' @importFrom splines ns
#' @export
#'
#'
get_index_comb <- function(fit, hurdle_fit = NULL, year = NULL, probs = c(0.025, 0.975), rescale = 1, pred_grid = NULL,  format = 'long', ...) {
  
  if  (!inherits(fit, c("sdmTMB", "glm", "survreg", "brmsfit"))) stop("This model class is not supported.")
  
  #_____________________________________________________________________________
  # Get all the necessary variables
  #_____________________________________________________________________________
  is_sdm <- inherits(fit, 'sdmTMB')
  
  # Name if first term
  if (is.null(year)) {
    year <- get_first_term(fit = fit)
  }
  # levels of fyear
  yrs <-  sort(unique(model.frame(fit)[[year]]))
  
  # number of years
  n <- length(yrs)
  
  #_____________________________________________________________________________
  #  Draw from model predictions for sdmTMB
  #_____________________________________________________________________________
  
  if(is_sdm) {
    # Spatio-temporal models 
    draws <- list()
    # Binomial Index 
    # come back to test this
    #  model <- if (length(spatiotemporal$formula)==2)  = 1
    predict_sim <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE, nsim = 1000, model = 1, type = "response")
    # generate draws of indices
    draws$Binomial <- sdmTMB::get_index_sims(predict_sim, return_sims = T) %>% # need to think about the area here
      rename(level = !!sym(year)) 
    
    #  Positive Index  
    predict_sim <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE, nsim = 1000, model = 2)
    # generate draws of indices
    draws$Positive <- sdmTMB::get_index_sims(predict_sim, return_sims = T) %>% # need to think about the area here
      rename(level = !!sym(year)) 
    
    
    
  } else {
    
    # create a list to hold binomial and positive model components if they exist
    mods <- list()
    mods$Positive <- if (is.null(fit$family$family) || !any(fit$family$family %in% c("bernoulli", "binomial"))) fit
    mods$Binomial <- if (!is.null(hurdle_fit)) {
      hurdle_fit
    } else if ((!is.null(fit$family$family) && any(fit$family$family %in% c("bernoulli", "binomial")))){
      fit
    }
    
    draws <- setNames(lapply(names(mods), function(component) {
      
      this_model <- mods[[component]]
      # model data
      raw_data <- this_model$data
      
      if (is.null(raw_data) && inherits(this_model, "survreg")) {
        raw_data <- eval(this_model$call$data)
      }
      
      # Forumula
      mod_terms <- all.vars(as.formula(formula(this_model)))
      
      # name of the response variable
      response_name <- insight::find_response(this_model)
      
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
      # Draw from model predictions for brms, GLM and survreg
      #_____________________________________________________________________________
      if (is.brmsfit(this_model)) {
        
        draws <- fitted(object = this_model, newdata = newdata, probs = c(probs[1], 0.5, probs[2]), re_formula = NA, scale = "response", summary = F)
        
        colnames(draws) <- yrs
        
        draws <- draws %>%
          as.data.frame() %>%
          mutate(.iteration = row_number()) %>%
          pivot_longer(cols = -.iteration,
                       names_to = 'level',
                       values_to = '.value') 
        
        
      }else if(inherits(fit, c("glm", "survreg"))){
        # GLMs and Survreg. 
        
        # Extract model coeffs
        cfs <- coefficients(this_model)
        
        # Extract covariance matrix. Survreg has extra param Log(scale), we subset vcov matrix to the size fo coefficients only 
        V <- vcov(this_model)[1:length(cfs), 1:length(cfs)]
        
        # Different process for binomial component
        if (!is.null(this_model$family$family) && any(this_model$family$family %in% c("bernoulli", "binomial"))){
          
          rows = c(1,which(substr(row.names(V),
                                  1,
                                  nchar(year))==year))
          
          resp_name <- names(model.frame(this_model))[1]
          int <- logit(aggregate(list(unstan=this_model$model[[resp_name]]),
                                 list(level=this_model$model[[year]]),
                                 mean)$unstan)[1]
          cfs[1] <- int
          
          bin <- mvtnorm::rmvnorm(1000,cfs,V)[,rows]
          bin <- cbind(inv_logit(bin[,1]),inv_logit(bin[,1] + bin[,2:length(yrs)]))
          
          draws <- bin
        } else {
          # Positive model component
          
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
          X <- model.matrix(delete.response(terms(this_model)), data = newdata)
          
          # Predicted response for new data
          draws <- exp(beta_draws %*% t(X))
          
        }
        
        # Rearrange df for index derivation
        colnames(draws) <- yrs
        
        draws <- draws %>%
          as.data.frame() %>%
          mutate(.iteration = row_number()) %>%
          pivot_longer(cols = -.iteration,
                       names_to = 'level',
                       values_to = '.value') 
        
      }
      return(draws)
    }), names(mods))
  }
  
  #________________________________________________________________________________
  # Bring it all together: Derive indices normalised with geometric mean
  #_________________________________________________________________________________
  # This step is common to all models
  
  # Derive Combined index if two models exist.
  if(length(draws)==2){
    draws$Combined <- draws$Positive
    draws$Combined$.value <- draws$Positive$.value * draws$Binomial$.value
  }
  
  # Do geometric mean transformation to draws, then summarise to find medians and quantiles
  
  indices <- setNames(lapply(names(draws), function(idx){
    index_stan <- draws[[idx]] %>%
      
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
    
    # Unstan indices
    if (idx != 'Combined'){
      predictor = ifelse(idx=='Binomial', 1, 2)
      this_model <- if (is_sdm) fit else mods[[idx]]
      indices <- get_unstandardised(fit = this_model, predictor = predictor)
      
      # merge unstan andf stan indices
      indices <- indices %>%
        left_join(index_stan,  by = 'level')
    } else { indices <- index_stan}
    return(indices)
  } ), names(draws))
  
  
  indices_wide <- bind_rows(indices, .id = "Index")
  
  if(format == 'wide'){
  return(indices_wide)
  }
  else{
    indices_long <- indices_wide %>%
      pivot_longer(
        cols = -c(level, Index),
        names_to = c("is_stan", "stat", "is_scaled"),
        names_pattern = "^(stan|unstan)(Lower|Upper)?_?(unscaled)?$",
        values_drop_na = TRUE
      ) %>%
      mutate(
        stat = ifelse(stat == "" , "median", stat),
        is_scaled = (is_scaled == ""),
        is_stan = (is_stan=='stan')
        # is_stan = factor(is_stan, levels = c("stan", "unstan"), 
        #                  labels = c("Standardised", "Unstandardised"))
        
      ) %>%
      # Pivot 'median', 'Lower', and 'Upper' back into their own columns
      pivot_wider(
        names_from = stat, 
        values_from = value
      )
    return(indices_long)
  }
}
