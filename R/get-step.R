#' Get the standardised indices with terms consecutively added
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
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @export
#' 
get_step <- function(fit, pred_grid = NULL, predictor = NULL) {
  
  if  (!inherits(fit, c("sdmTMB", "glm", "survreg", "brmsfit"))) stop("This model class is not supported.")
  
  is_sdm <- inherits(fit, 'sdmTMB')
  
  if (is_sdm){
    
    if  (is.null(predictor)) stop("Argument 'predictor' is missing. Please specify 1 for the first part or 2 for the hurdle part.")
    if  (is.null(pred_grid)) stop("Argument 'pred_grid' is missing. Please pass prediction grid including environment covariate values for the full geograhpical area.")
    newFormula <- (fit$formula)[[predictor]]
    sptp_on <- fit$spatiotemporal[predictor]!="off"
    
  } else {
    # brms, GLM and survreg
    newFormula <- formula(fit)
    sptp_on <- FALSE
    
  }
  
  # extract terms from this formula
  terms_labels <- get_terms(fit, predictor = predictor)
  
  # initiate effects list
  effects <- list()
  
  # initiate summary table
  step_summary <- list()
  
  # Create models with terms successively added
  for(termCount in 0:(length(terms_labels)+is_sdm + sptp_on)){
    
    if(termCount>0){
      
      term <- terms_labels[termCount]
      
      # Update both formula and model
      
      newFormula <- update.formula(newFormula,
                                   formula(paste("~",
                                                 paste(paste(terms_labels[1:min(termCount,length(terms_labels)) ],
                                                             collapse='+')))))
      
      # Prepare arguments for the re-fit
      fit_args <- list(object = fit, formula. = newFormula)
      
      # Conditionally add spatiotemoral terms to model call
      if(is_sdm){
        # Turn spatial and/or spatio-temporal component on
        # spatial on only if this iteration is after all terms were added
        fit_args$spatial <- ifelse(termCount<=length(terms_labels), 'off', fit$spatial[predictor])
        # spatiotemporal 'on' only if this is the last iteration AND original fit has it 'on'
        fit_args$spatiotemporal <- ifelse(termCount!=length(terms_labels)+2, 'off', fit$spatiotemporal[predictor])
      }
      
      # refit the model
      fit_reduced <- do.call(update, fit_args)
      
      # Get index for this model
      idx_reduced <- get_index_comb (fit_reduced, pred_grid = pred_grid, predictor = predictor)
      
      # print(summary(fit_reduced))
      # Generate the right hand side of formula as name for index
      idx_name <- case_when(
        termCount == 1                           ~ term,
        termCount == length(terms_labels) + 2    ~ "+ spatiotemporal",
        termCount == length(terms_labels) + 1    ~ "+ spatial",
        TRUE                                     ~ paste("+", term)
      )
      
      # Store column of indices
      effects[[idx_name]] <- idx_reduced %>%
        filter(is_stan, is_scaled)%>% select(Index, level, median, Lower, Upper)
      
    } else {
      term = 'intercept'
      fit_reduced = update(fit,.~1)
    }
    
    # TO DO calculate summary statistics here
    
    logLik <- switch(class(fit_reduced)[1],
                     brmsfit = mean(rowSums(log_lik(fit_reduced))),
                     logLik(fit_reduced))
    
    aic <- switch(class(fit_reduced)[1],
                  brmsfit = NA,
                  AIC(fit_reduced))
    
    r2Dev <- switch(class(fit_reduced)[1],
                    glm = (fit_reduced$null.deviance-fit_reduced$deviance)/fit_reduced$null.deviance*100,
                    NA)
    
    step_summary[[termCount+1]] <- data.frame(term = term,
                                              df = extractAIC(fit_reduced)[1],
                                              logLik = logLik,
                                              AIC = aic,
                                              r2Dev = r2Dev)
    
  }
  
  # combine summary stats:
  step_summary <- do.call(rbind, step_summary)  
  
  
  if (inherits(fit, 'survreg')) {
    
    resp <- as.character(formula(fit)[2])
    resp <- gsub("^Surv\\(|\\)$", "", resp)
    saturated_model <- update(fit, formula. = paste('. ~ ',resp))
    residDev <- 2*(saturated_model$loglik[2] - step_summary$logLik)
    step_summary$r2Dev <- (residDev[1] - residDev)/residDev[1]*100
    
  }
  
  step_summary %<>% 
    mutate(
      r2Dev_delta = r2Dev - lag(r2Dev, default =  0), # defaults are there to handle first row. It is the value of lag at 'time zero'
      df = df - lag(df, default = 0),
      Included = "*"
    ) %>%
    select(-logLik)
  
  #  Combine all the effect columns into one wide data frame
  all_idx <- bind_rows(effects, .id = "Model") %>%
    mutate(Model = factor(Model, levels = unique(Model)))
  
  #  Final bind indices from last iteration with unstan and CI + all the steps 
  # This was ghoti process, but it does not seem to be used anywhere
  # indices <- cbind(idx_reduced, all_idx)  
  
  return(list(step_indices = all_idx, step_summary = step_summary))
}

