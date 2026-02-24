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
#' @importFrom stats fitted
#' @importFrom brms is.brmsfit
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
      idx_reduced <- get_index (fit_reduced,  pred_grid = pred_grid, predictor = predictor)
      # print(summary(fit_reduced))
      # Generate the right hand side of formula as name for index
      idx_name <- case_when(
        termCount == 1                           ~ term,
        termCount == length(terms_labels) + 2    ~ "+ spatiotemporal",
        termCount == length(terms_labels) + 1    ~ "+ spatial",
        TRUE                                     ~ paste("+", term)
      )
      
      # Store column of indices
      effects[[idx_name]] <- idx_reduced$stan
      
    } else {
      term = 'intercept'
      fit_reduced = update(fit,.~1)
    }
    
    # TO DO calculate summary statistics here
    
  }
  #  Combine all the effect columns into one wide data frame
  all_idx <- do.call(cbind, effects)
  
  #  Final bind indices from last iteration with unstan and CI + all the steps 
  indices <- cbind(idx_reduced, all_idx)  
  
  return(indices)
}


