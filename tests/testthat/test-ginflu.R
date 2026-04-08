# Tests for the Ginflu package
#
# Primary goal: verify that data structures returned by Ginflu functions are
# consistent across supported model types (glm gaussian, poisson, binomial,
# survreg) and that output formats are correct.

library(Ginflu)

# ---------------------------------------------------------------------------
# Shared test data
# ---------------------------------------------------------------------------

set.seed(42)
n <- 60
test_df <- data.frame(
  year   = factor(rep(2001:2006, each = 10)),
  month  = factor(rep(c("A", "B"), times = 30)),
  effort = runif(n, 0.5, 2),
  catch  = exp(
    0.3 * rep(1:6, each = 10) +
    0.2 * rep(c(0, 1), 30) +
    rnorm(n, 0, 0.3)
  ),
  # Use a generous lambda to avoid zeros in count (log(0) produces NA in
  # get_unstandardised when the positive component is used).
  count   = rpois(n, lambda = 5 + rep(1:6, each = 10)),
  present = rbinom(n, 1, prob = plogis(0.5 * rep(1:6, each = 10) / 3 - 0.5))
)

# Bake the data frame into the model call via bquote() so that update() inside
# get_step() does not need to search for 'test_df' by name when running inside
# the testthat environment.
glm_gauss   <- eval(bquote(glm(log(catch) ~ year + month + effort, data = .(test_df), family = gaussian())))
glm_poisson <- eval(bquote(glm(count ~ year + month,               data = .(test_df), family = poisson())))
glm_binom   <- eval(bquote(glm(present ~ year + month,             data = .(test_df), family = binomial())))

n_years <- length(levels(test_df$year))


# ===========================================================================
# Helper utilities
# ===========================================================================

test_that("inv_logit returns values in (0, 1)", {
  x <- c(-10, -1, 0, 1, 10)
  out <- inv_logit(x)
  expect_true(all(out > 0 & out < 1))
  expect_equal(inv_logit(0), 0.5)
})

test_that("gmean returns geometric mean", {
  expect_equal(gmean(c(1, 1, 1)), 1)
  expect_equal(gmean(c(2, 8)), 4)       # geometric mean of 2 and 8 is 4
  expect_equal(gmean(c(1, 4)), 2)
})

test_that("mean_or_mode returns mean for numeric", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(mean_or_mode(x), mean(x))
})

test_that("mean_or_mode returns a factor for character/factor input", {
  x <- c("a", "a", "b", "c")
  out <- mean_or_mode(x)
  expect_true(is.factor(out))
  expect_equal(as.character(out), "a")   # mode is "a"
})


# ===========================================================================
# get_first_term
# ===========================================================================

test_that("get_first_term returns 'year' for glm_gauss", {
  expect_equal(get_first_term(glm_gauss), "year")
})

test_that("get_first_term returns 'year' for glm_poisson", {
  expect_equal(get_first_term(glm_poisson), "year")
})

test_that("get_first_term returns 'year' for glm_binom", {
  expect_equal(get_first_term(glm_binom), "year")
})


# ===========================================================================
# get_terms
# ===========================================================================

test_that("get_terms returns a character vector for glm_gauss", {
  out <- get_terms(glm_gauss)
  expect_true(is.character(out))
  expect_true("year" %in% out)
  expect_true("month" %in% out)
  expect_true("effort" %in% out)
})

test_that("get_terms returns a character vector for glm_poisson", {
  out <- get_terms(glm_poisson)
  expect_true(is.character(out))
  expect_true("year" %in% out)
  expect_true("month" %in% out)
})

test_that("get_terms returns same length vector for glm_gauss and glm_poisson (3 vs 2 terms)", {
  expect_equal(length(get_terms(glm_gauss)), 3)
  expect_equal(length(get_terms(glm_poisson)), 2)
})

test_that("get_terms structure is consistent across glm types", {
  for (fit in list(glm_gauss, glm_poisson, glm_binom)) {
    out <- get_terms(fit)
    expect_true(is.character(out))
    expect_true(length(out) >= 1)
  }
})


# ===========================================================================
# get_preds
# ===========================================================================

test_that("get_preds returns a list with required components for glm_gauss", {
  out <- get_preds(glm_gauss)
  expect_true(is.list(out))
  expect_named(out, c("preds", "V", "X_centered", "assign", "terms", "year"),
               ignore.order = TRUE)
})

test_that("get_preds$preds is a data.frame with fit.* and se.fit.* columns", {
  out <- get_preds(glm_gauss)
  preds_df <- out$preds
  expect_true(is.data.frame(preds_df))

  terms <- get_terms(glm_gauss)
  for (trm in terms) {
    expect_true(paste0("fit.", trm) %in% names(preds_df),
                info = paste("Missing fit.", trm))
    expect_true(paste0("se.fit.", trm) %in% names(preds_df),
                info = paste("Missing se.fit.", trm))
  }
})

test_that("get_preds$V is a square matrix for glm_gauss", {
  out <- get_preds(glm_gauss)
  V <- out$V
  expect_true(is.matrix(V))
  expect_equal(nrow(V), ncol(V))
})

test_that("get_preds$terms matches get_terms for glm_gauss", {
  out <- get_preds(glm_gauss)
  expect_equal(out$terms, get_terms(glm_gauss))
})

test_that("get_preds$year is 'year' for all glm types", {
  for (fit in list(glm_gauss, glm_poisson, glm_binom)) {
    out <- get_preds(fit)
    expect_equal(out$year, "year")
  }
})

test_that("get_preds returns consistent list structure across glm types", {
  required_names <- c("preds", "V", "X_centered", "assign", "terms", "year")
  for (fit in list(glm_gauss, glm_poisson, glm_binom)) {
    out <- get_preds(fit)
    expect_true(all(required_names %in% names(out)),
                info = paste("Missing components for", class(fit)[1]))
  }
})

test_that("get_preds$preds has same number of rows as data", {
  out <- get_preds(glm_gauss)
  expect_equal(nrow(out$preds), nrow(test_df))
})


# ===========================================================================
# get_unstandardised
# ===========================================================================

test_that("get_unstandardised returns a data.frame for glm_gauss", {
  out <- get_unstandardised(glm_gauss)
  expect_true(is.data.frame(out))
})

test_that("get_unstandardised has 'level' and 'unstan' columns for positive models", {
  for (fit in list(glm_gauss, glm_poisson)) {
    out <- get_unstandardised(fit)
    expect_true("level" %in% names(out), info = class(fit)[1])
    expect_true("unstan" %in% names(out), info = class(fit)[1])
  }
})

test_that("get_unstandardised returns one row per year level", {
  for (fit in list(glm_gauss, glm_poisson)) {
    out <- get_unstandardised(fit)
    expect_equal(nrow(out), n_years, info = class(fit)[1])
  }
})

test_that("get_unstandardised unstan values are positive", {
  for (fit in list(glm_gauss, glm_poisson)) {
    out <- get_unstandardised(fit)
    expect_true(all(out$unstan > 0), info = class(fit)[1])
  }
})

test_that("get_unstandardised for binomial has 'unstan' and 'unstan_unscaled' columns", {
  out <- get_unstandardised(glm_binom)
  expect_true("level" %in% names(out))
  expect_true("unstan" %in% names(out))
  expect_true("unstan_unscaled" %in% names(out))
})

test_that("get_unstandardised for binomial returns one row per year level", {
  out <- get_unstandardised(glm_binom)
  expect_equal(nrow(out), n_years)
})

test_that("get_unstandardised unstan values are positive for binomial", {
  out <- get_unstandardised(glm_binom)
  expect_true(all(out$unstan > 0))
})


# ===========================================================================
# get_index — structure and column consistency across model types
# ===========================================================================

# Compute indices once (expensive) and reuse across tests
idx_gauss   <- get_index(glm_gauss)
idx_poisson <- get_index(glm_poisson)
idx_binom   <- get_index(glm_binom)

test_that("get_index returns a data.frame for all glm types", {
  for (idx in list(idx_gauss, idx_poisson, idx_binom)) {
    expect_true(is.data.frame(idx))
  }
})

test_that("get_index returns one row per year level for all glm types", {
  for (idx in list(idx_gauss, idx_poisson, idx_binom)) {
    expect_equal(nrow(idx), n_years)
  }
})

test_that("get_index has required standardised-index columns for positive models", {
  required <- c("level", "unstan", "stan", "stanLower", "stanUpper",
                "stan_unscaled", "stanLower_unscaled", "stanUpper_unscaled")
  for (idx in list(idx_gauss, idx_poisson)) {
    expect_true(all(required %in% names(idx)),
                info = paste("Missing:", paste(setdiff(required, names(idx)), collapse = ", ")))
  }
})

test_that("get_index has required columns for binomial model", {
  required <- c("level", "unstan", "stan", "stanLower", "stanUpper",
                "stan_unscaled", "stanLower_unscaled", "stanUpper_unscaled")
  expect_true(all(required %in% names(idx_binom)),
              info = paste("Missing:", paste(setdiff(required, names(idx_binom)), collapse = ", ")))
})

test_that("get_index column names are the same for gaussian and poisson", {
  expect_equal(sort(names(idx_gauss)), sort(names(idx_poisson)))
})

test_that("get_index stan values are positive for all glm types", {
  for (idx in list(idx_gauss, idx_poisson, idx_binom)) {
    expect_true(all(idx$stan > 0))
    expect_true(all(idx$stanLower > 0))
    expect_true(all(idx$stanUpper > 0))
  }
})

test_that("get_index stanLower <= stan <= stanUpper for all glm types", {
  for (idx in list(idx_gauss, idx_poisson, idx_binom)) {
    expect_true(all(idx$stanLower <= idx$stan + 1e-8))
    expect_true(all(idx$stan <= idx$stanUpper + 1e-8))
  }
})

test_that("get_index stan has geometric mean approximately equal to 1", {
  # get_index uses median(rel_idx) per year; the geometric mean of those medians
  # is close but not exactly 1 due to Jensen's inequality, so tolerance is
  # relaxed to 5 %.
  for (idx in list(idx_gauss, idx_poisson, idx_binom)) {
    gm <- exp(mean(log(idx$stan)))
    expect_equal(gm, 1, tolerance = 0.05)
  }
})

test_that("get_index level column matches year factor levels", {
  for (idx in list(idx_gauss, idx_poisson, idx_binom)) {
    expect_true(all(as.character(idx$level) %in% as.character(levels(test_df$year))))
  }
})


# ===========================================================================
# get_step — structure across model types
# ===========================================================================

step_gauss   <- get_step(glm_gauss)
step_poisson <- get_step(glm_poisson)

test_that("get_step returns a named list with 'step_indices' and 'step_summary'", {
  for (out in list(step_gauss, step_poisson)) {
    expect_true(is.list(out))
    expect_named(out, c("step_indices", "step_summary"), ignore.order = TRUE)
  }
})

test_that("get_step$step_summary is a data.frame with expected columns", {
  expected_cols <- c("term", "df", "AIC", "r2Dev", "r2Dev_delta", "Included")
  for (out in list(step_gauss, step_poisson)) {
    ss <- out$step_summary
    expect_true(is.data.frame(ss))
    expect_true(all(expected_cols %in% names(ss)),
                info = paste("Missing:", paste(setdiff(expected_cols, names(ss)), collapse = ", ")))
  }
})

test_that("get_step$step_summary has one row per term plus intercept", {
  # terms + intercept row
  expect_equal(nrow(step_gauss$step_summary),   length(get_terms(glm_gauss)) + 1)
  expect_equal(nrow(step_poisson$step_summary), length(get_terms(glm_poisson)) + 1)
})

test_that("get_step$step_summary first row is the intercept", {
  for (out in list(step_gauss, step_poisson)) {
    expect_equal(out$step_summary$term[1], "intercept")
  }
})

test_that("get_step$step_indices is a data.frame", {
  for (out in list(step_gauss, step_poisson)) {
    expect_true(is.data.frame(out$step_indices))
  }
})

test_that("get_step$step_indices contains all columns from get_index", {
  index_cols <- names(idx_gauss)
  si <- step_gauss$step_indices
  expect_true(all(index_cols %in% names(si)),
              info = paste("Missing:", paste(setdiff(index_cols, names(si)), collapse = ", ")))
})

test_that("get_step$step_indices has one row per year level", {
  for (out in list(step_gauss, step_poisson)) {
    expect_equal(nrow(out$step_indices), n_years)
  }
})

test_that("get_step$step_indices has step columns for each non-intercept term", {
  terms_g <- get_terms(glm_gauss)
  si <- step_gauss$step_indices
  # Step columns are the terms themselves (first term) or '+ term'
  step_cols <- c(terms_g[1], paste("+", terms_g[-1]))
  expect_true(all(step_cols %in% names(si)),
              info = paste("Missing step columns:",
                           paste(setdiff(step_cols, names(si)), collapse = ", ")))
})

test_that("get_step column structure is consistent across glm_gauss and glm_poisson", {
  # Both should have the same non-step base columns
  base_cols <- c("level", "unstan", "stan", "stanLower", "stanUpper")
  for (out in list(step_gauss, step_poisson)) {
    si <- out$step_indices
    expect_true(all(base_cols %in% names(si)))
  }
})

test_that("get_step step_summary r2Dev is non-decreasing (adding terms improves fit)", {
  # cumulative r2Dev should be non-decreasing (each term adds explanatory power)
  for (out in list(step_gauss, step_poisson)) {
    r2_vals <- out$step_summary$r2Dev
    r2_cumsum <- cumsum(replace(out$step_summary$r2Dev_delta,
                                is.na(out$step_summary$r2Dev_delta), 0))
    expect_true(all(diff(r2_cumsum) >= -1e-6))
  }
})


# ===========================================================================
# Error handling
# ===========================================================================

test_that("get_unstandardised errors for unsupported model class", {
  lm_fit <- lm(log(catch) ~ year + month, data = test_df)
  expect_error(get_unstandardised(lm_fit))
})

test_that("get_index errors for unsupported model class", {
  lm_fit <- lm(log(catch) ~ year + month, data = test_df)
  expect_error(get_index(lm_fit))
})

test_that("get_step errors for unsupported model class", {
  lm_fit <- lm(log(catch) ~ year + month, data = test_df)
  expect_error(get_step(lm_fit))
})
