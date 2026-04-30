#' Plot Coefficient Distribution and Influence (CDI)
#'
#' @description
#' Generates a multi-pane layout visualising the coefficients, data distribution,
#' and relative influence of a model predictor.
#'
#' @param fit A fitted model object (e.g., from `glm` or `gam`).
#' @param year Character string. The name of the temporal column (e.g., "fyear")
#'   to be used as the focus variable on the y-axis of the distribution plot.
#' @import ggplot2
#' @import dplyr
#' @importFrom patchwork plot_layout
#' @importFrom stats as.formula coef
#' @importFrom cowplot background_grid
#' @export

plot_cdi <- function(preds_list, compare_preds_list = NULL) {
  # comparison switch

  compareOn <- !is.null(compare_preds_list)

  if (compareOn) {
    compare_preds_df <- compare_preds_list$preds
  }

  # Extract components from the preds list

  preds <- preds_list$preds
  V <- preds_list$V
  X_centered <- preds_list$X_centered
  assigns <- preds_list$assign
  all_model_terms <- preds_list$terms
  year <- preds_list$year

  # exclude year from terms
  terms_labels <- setdiff(all_model_terms, year)

  result <- setNames(
    lapply(terms_labels, function(term_label) {
      # Define col names for fitted terms and ses
      fit_colname <- sym(paste0('fit.', term_label))
      se_colname <- sym(paste0('se.fit.', term_label))

      # Identify model matrix columns for this specific term
      match_idx <- which(all_model_terms == term_label)
      cols <- which(assigns == match_idx)

      # Extract the sub model matrix for this term
      Xi_centered <- X_centered[, cols, drop = FALSE]
      Vi <- V[cols, cols, drop = FALSE]

      # Define levels of term on which coefficient and distribution plots will be based
      # This is done by search for each column name in the data as a whole word in the
      # each term. This allows for matching of terms like 'poly(log(var),3)' with 'var'

      term_stripped <- all.vars(as.formula(paste("~", term_label)))
      raw_values <- preds[[term_stripped]]

      # Numeric terms are cut into factors
      if (
        is.numeric(raw_values) &&
          (!all(raw_values %% 1 == 0) | length(unique(raw_values)) > 10)
      ) {
        breaks <- pretty(raw_values, 30)
        step <- breaks[2] - breaks[1]
        breaks <- breaks[breaks <= max(raw_values)]

        if (max(breaks) < max(raw_values)) {
          breaks <- c(breaks, max(breaks) + step)
        }
        labels <- breaks[-length(breaks)] + (step / 2)
        levels <- cut(
          raw_values,
          breaks = breaks,
          labels = labels,
          include.lowest = TRUE
        )
      } else {
        levels <- factor(raw_values)
      }

      coeffs <- preds %>%
        mutate(row_id = row_number()) %>%
        group_by(term = !!levels) %>%
        summarise(
          coef = mean(!!fit_colname),
          se = {
            # Calculate the average design matrix row for this bin
            X_bin_avg <- colMeans(Xi_centered[row_id, , drop = FALSE])
            # Project the covariance matrix: sqrt(avg_row %*% V %*% avg_row_T)
            sqrt(t(X_bin_avg) %*% Vi %*% X_bin_avg)
          },
          .groups = "drop"
        ) %>%
        mutate(
          lower = coef - (1.96 * se),
          upper = coef + (1.96 * se)
        )

      # Reorder levels according to coefficients for factors
      if (is.factor(raw_values) && !(grepl('month', term_label))) {
        coeffs <- coeffs %>%
          arrange(coef)

        coeffs$term <- factor(coeffs$term, levels = coeffs$term, ordered = TRUE)

        levels <- factor(levels, levels = coeffs$term, ordered = TRUE)
      }

      if (
        compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
      ) {
        # Extract components from the comparison list
        comp_raw_values <- compare_preds_df[[term_stripped]]
        comp_preds <- compare_preds_list$preds
        comp_V <- compare_preds_list$V
        comp_X_centered <- compare_preds_list$X_centered
        comp_assign <- compare_preds_list$assign
        comp_terms <- compare_preds_list$terms

        # Numeric terms are cut into factors
        if (
          is.numeric(comp_raw_values) &&
            (!all(comp_raw_values %% 1 == 0) |
              length(unique(comp_raw_values)) > 10)
        ) {
          comp_levels <- cut(
            comp_raw_values,
            breaks = breaks,
            labels = labels,
            include.lowest = TRUE
          )
        } else {
          comp_levels <- factor(comp_raw_values)
        }

        # Identify columns in the comparison model matrix
        comp_match_idx <- which(comp_terms == term_label)
        comp_cols <- which(comp_assign == comp_match_idx)
        # Extract the sub model matrix for this term
        comp_Xi_centered <- comp_X_centered[, comp_cols, drop = FALSE]
        comp_Vi <- comp_V[comp_cols, comp_cols, drop = FALSE]

        comp_coeffs <- compare_preds_df %>%
          mutate(row_id = row_number()) %>%
          group_by(term = !!comp_levels) %>%
          summarise(
            coef = mean(!!fit_colname),
            se = {
              X_bin_avg <- colMeans(comp_Xi_centered[row_id, , drop = FALSE])
              as.numeric(sqrt(t(X_bin_avg) %*% comp_Vi %*% X_bin_avg))
            },
            .groups = "drop"
          ) %>%
          mutate(
            term = factor(term, levels = levels(coeffs$term), ordered = TRUE),
            lower = coef - (1.96 * se),
            upper = coef + (1.96 * se)
          )
      }

      # (1) Coefficients  Plot
      #_____________________________________________________________________________

      # Colours when comparison, black otherwise?
      col1 <- if (compareOn) 'dodgerblue1' else 'black'

      # Create plot theme depending on the term:

      if ((length(levels(levels)) > 12)) {
        # Vertical labels and tighter margins
        dynamic_theme <- theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.x.top = element_text(vjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 15))
        )
      } else {
        # Horizontal labels
        dynamic_theme <- theme(
          axis.text.x = element_text(angle = 0),
          axis.title.x = element_text(margin = margin(t = 10))
        )
      }

      # Create theme elements used for each of the 3 plots
      common_theme <- theme(
        axis.title.x = element_text(vjust = 0),
        panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
        axis.line.x = element_line(colour = "black", linewidth = 0.6),
        axis.line.y = element_line(colour = "black", linewidth = 0.6)
      )

      # No x-labels for vessel_key
      x_labels <- if (grepl('vessel_key', term_label, ignore.case = TRUE)) {
        NULL
        # }  else if(length(levels(levels)) > 15) {
        #   # Downsample labels for continuous vars otherwise they overlap, the logic is to keep every n-th label where n = original no of levels/10
        #   function(x) {
        #     x[seq_along(x) %% max(1, round(length(x) / 10)) != 1] <- ""
        #     return(x)
        #   }
      } else {
        waiver()
      }

      p1 <- ggplot(coeffs, aes(x = term, y = exp(coef))) +
        geom_point(shape = 2, size = 2, color = col1) +
        geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), width = 0) +
        # Bottom Cap
        geom_segment(aes(
          x = as.numeric(term) - 0.05,
          xend = as.numeric(term) + 0.05,
          y = exp(lower),
          yend = exp(lower)
        )) +
        # Top Cap
        geom_segment(aes(
          x = as.numeric(term) - 0.05,
          xend = as.numeric(term) + 0.05,
          y = exp(upper),
          yend = exp(upper)
        )) +
        geom_hline(yintercept = 1, linetype = "dashed") +

        # Conditional block to display coeffs for a model to be compared with (e.g. last year)
        {
          if (
            compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
          ) {
            list(
              geom_point(
                data = comp_coeffs,
                aes(x = term, y = exp(coef)),
                shape = 17,
                size = 2,
                color = '#E41A1CCC'
              ),
              geom_errorbar(
                data = comp_coeffs,
                aes(ymin = exp(lower), ymax = exp(upper)),
                color = '#E41A1CCC',
                width = 0
              )
            )
          }
        } +

        # Trying to place enough labels on the y-scale
        # Solution 1:
        scale_y_log10(breaks = scales::breaks_extended(n = 9)) +

        # Solution 2:
        # {
        #   if (max(exp(coeffs$coef)) - min(exp(coeffs$coef)) > 1) {
        #     scale_y_log10(breaks = scales::breaks_log(n = 8, base = 1.2),
        #                   labels = scales::label_number(accuracy = 0.01))
        #   } else {
        #     scale_y_log10(n.breaks = 8)
        #   }
        # } +
        scale_x_discrete(drop = FALSE, labels = x_labels, position = 'top') +
        theme_bw() +
        background_grid(major = "x", minor = "none") +
        labs(y = 'Effect', x = NULL) +
        dynamic_theme +
        common_theme

      # (2) Distribution Plot
      #_____________________________________________________________________________

      distrs <- preds %>%
        group_by(focus = .[[year]], term = levels) %>%
        summarise(count = n(), .groups = "drop_last") %>%

        mutate(
          total = sum(count),
          prop = count / total
        ) %>%

        ungroup()

      if (
        compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
      ) {
        comp_distrs <- comp_preds %>%
          group_by(focus = .[[year]], term = comp_levels) %>%
          summarise(count = n(), .groups = "drop_last") %>%

          mutate(
            total = sum(count),
            prop = count / total
          ) %>%

          ungroup()
      }

      p2 <- ggplot(distrs, aes(x = term, y = focus, size = sqrt(prop) * 20)) +
        geom_point(pch = 1, col = col1) +

        # Conditional block to display dists for a model to be compared with (e.g. last year)
        {
          if (
            compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
          ) {
            geom_point(
              data = comp_distrs,
              pch = 1,
              color = '#E41A1CCC',
              alpha = 0.6
            )
          }
        } +

        scale_size_identity() +
        scale_x_discrete(drop = FALSE, labels = x_labels) +
        theme_bw() +
        background_grid(major = "xy", minor = "none") +
        labs(x = term_stripped, y = year) +
        dynamic_theme +
        common_theme

      # (3) Influence
      #______________________________________________________________________________

      # This is formula 3b from Bentley et al. 2012
      infl <- preds %>%
        group_by(level = !!sym(year)) %>%
        summarise(Influence = mean(!!sym(fit_colname)))

      # Mean absolute magnitude of the term's impact, converted to annual % change
      overall <- exp(mean(abs(infl$Influence))) - 1

      # Annual growth % of this term's contribution
      n <- nrow(infl)
      trend <- exp(cov(1:n, infl$Influence) / var(1:n)) - 1

      infl <- infl %>%
        mutate(
          Influence = exp(Influence),
          Term = term_label,
          overall = overall,
          trend = trend
        ) %>%
        relocate(Term, .after = level)

      if (
        compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
      ) {
        comp_infl <- compare_preds_df %>%
          group_by(level = !!sym(year)) %>%
          summarise(Influence = exp(mean(!!sym(fit_colname))))
        comp_x_min <- min(comp_infl$Influence)
        comp_x_max <- max(comp_infl$Influence)
      }

      p3 <- ggplot(infl, aes(x = Influence, y = level)) +

        # Conditional block to display coeffs for a model to be compared with (e.g. last year)
        # In the background: The two rectangles (if compareOn is TRUE)
        {
          if (
            compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
          ) {
            list(
              # Outer shaded rectangle (+/- 0.05 padding)
              annotate(
                "rect",
                xmin = comp_x_min - 0.05,
                xmax = comp_x_max + 0.05,
                ymin = -Inf,
                ymax = Inf,
                fill = "grey90",
                alpha = 0.5
              ),
              # Inner shaded rectangle (exact range)
              annotate(
                "rect",
                xmin = comp_x_min,
                xmax = comp_x_max,
                ymin = -Inf,
                ymax = Inf,
                fill = "grey80",
                alpha = 0.5
              )
            )
          }
        } +

        # In the middle layer: Black dots and lines for current (or only) model
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_line(group = 1, col = col1) +
        geom_point(size = 3, pch = 16, col = col1) +

        # In the foreground: Brown dots and lines for model to compare with
        {
          if (
            compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
          ) {
            list(
              geom_line(
                data = comp_infl,
                group = 1,
                linetype = 'dashed',
                color = '#E41A1CCC',
                alpha = 0.6
              ),
              geom_point(
                data = comp_infl,
                shape = 16,
                size = 3,
                color = '#E41A1CCC',
                alpha = 0.6
              )
            )
          }
        } +

        scale_x_continuous(
          limits = function(x) {
            c(min(pretty(x), 0.8), max(pretty(x), 1.2))
          },
          labels = scales::label_number(accuracy = 0.01)
        ) +
        scale_y_discrete(position = "right") +
        labs(x = "Influence", y = NULL) +
        theme_bw() +
        background_grid(major = "y", minor = "none") +
        common_theme

      # (4) Blank or legend for top right corner
      #______________________________________________________________________________
      pv <- ggplot() +
        geom_point(aes(x = 0, y = 0), alpha = 0) +

        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_line(colour = "black", linewidth = 0.6),
          axis.line.y = element_line(colour = "black", linewidth = 0.6)
        )

      # Legend for comparison plot

      if (
        compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
      ) {
        labels = c(
          "Current",
          paste(
            "Previous update to",
            max(as.numeric(as.character(comp_infl$level)))
          ),
          "Range of previous update",
          "Range of \n previous update +5%"
        )

        pv <- ggplot() +
          # Row 1: Current
          annotate(
            "segment",
            x = 5,
            xend = 15,
            y = 50,
            yend = 50,
            color = "dodgerblue1",
            linetype = "solid"
          ) +
          annotate(
            "point",
            x = 10,
            y = 50,
            shape = 16,
            size = 2,
            color = "dodgerblue1"
          ) +
          annotate(
            "point",
            x = 3,
            y = 50,
            shape = 2,
            size = 2,
            color = "dodgerblue1"
          ) +

          # Row 2: Previous
          annotate(
            "segment",
            x = 5,
            xend = 15,
            y = 40,
            yend = 40,
            color = "#E41A1CCC",
            linetype = "dotted"
          ) +
          annotate(
            "point",
            x = 10,
            y = 40,
            shape = 16,
            size = 2,
            color = "#E41A1CCC"
          ) +
          annotate(
            "point",
            x = 3,
            y = 40,
            shape = 17,
            size = 2,
            color = "#E41A1CCC"
          ) +

          # Row 3 & 4: Ranges
          annotate(
            "point",
            x = 10,
            y = 30,
            shape = 15,
            size = 6,
            color = "grey80",
            alpha = 0.5
          ) +
          annotate(
            "point",
            x = 10,
            y = 20,
            shape = 15,
            size = 6,
            color = "grey90",
            alpha = 0.5
          ) +

          # Text Labels
          annotate(
            "text",
            x = 16,
            y = seq(from = 50, to = 20, by = -10),
            label = labels,
            hjust = 0,
            size = 4
          ) +

          coord_cartesian(xlim = c(0, 60), ylim = c(0, 60), expand = FALSE) +

          theme_minimal() +
          theme(
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line.x = element_line(colour = "black", linewidth = 0.6),
            axis.line.y = element_line(colour = "black", linewidth = 0.6)
          )
      }

      combined_plot <- p1 +
        pv +
        p2 +
        p3 +
        plot_layout(nrow = 2, ncol = 2, heights = c(1, 2), widths = c(2, 1)) &
        theme(
          # axis.title = element_text(size = 14),
          #  axis.text = element_text(size = 14),
          plot.margin = margin(0, 0, 0, 0)
        )

      # print(combined_plot)
    }),
    terms_labels
  )

  return(result)
}


#' Plot Coefficient Distribution and Influence (CDI) and calculate comparison traffic light metrics
#'
#' @description
#' Generates a multi-pane layout visualising the coefficients, data distribution,
#' and relative influence of a model predictor.
#'
#' @param fit A fitted model object (e.g., from `glm` or `gam`).
#' @param year Character string. The name of the temporal column (e.g., "fyear")
#'   to be used as the focus variable on the y-axis of the distribution plot.
#' @import ggplot2
#' @import dplyr
#' @importFrom patchwork plot_layout
#' @importFrom stats as.formula coef
#' @importFrom cowplot background_grid
#' @export
#'
cdi_plot_with_indicators <- function(preds_list, compare_preds_list = NULL) {
  
  # ------------------- Prepare the data ----------------------
  # comparison switch
  compareOn <- !is.null(compare_preds_list)

  if (compareOn) {
    compare_preds_df <- compare_preds_list$preds
  }

  indicators <- NULL # initialise indicators dataframe populateted if compareOn

  # Extract components from the preds list

  preds <- preds_list$preds
  V <- preds_list$V
  X_centered <- preds_list$X_centered
  assigns <- preds_list$assign
  all_model_terms <- preds_list$terms
  year <- preds_list$year

  # exclude year from terms
  terms_labels <- setdiff(all_model_terms, year)

  # ---------- Loop through terms to generate plots ------------

  result <- setNames(
    lapply(terms_labels, function(term_label) {

      # ---------- Organise data for this term ------------------
      # Define col names for fitted terms and standard errors
      fit_colname <- sym(paste0('fit.', term_label))
      se_colname <- sym(paste0('se.fit.', term_label))

      # Identify model matrix columns for this specific term
      match_idx <- which(all_model_terms == term_label)
      cols <- which(assigns == match_idx)

      # Extract the sub model matrix for this term
      Xi_centered <- X_centered[, cols, drop = FALSE]
      Vi <- V[cols, cols, drop = FALSE]

      # Define levels of term on which coefficient and distribution plots will be based
      # This is done by search for each column name in the data as a whole word in the
      # each term. This allows for matching of terms like 'poly(log(var),3)' with 'var'

      term_stripped <- all.vars(as.formula(paste("~", term_label)))
      if(length(term_stripped)==2) {
        raw_values <- paste(preds[[ term_stripped[1] ]], 
                   preds[[ term_stripped[2] ]], 
                   sep = ":")

      } else {
      raw_values <- preds[[term_stripped]]
      }

      # Numeric terms are cut into factors
      if (
        is.numeric(raw_values) &&
          (!all(raw_values %% 1 == 0) | length(unique(raw_values)) > 10)
      ) {
        breaks <- pretty(raw_values, 30)
        step <- breaks[2] - breaks[1]
        breaks <- breaks[breaks <= max(raw_values)]

        if (max(breaks) < max(raw_values)) {
          breaks <- c(breaks, max(breaks) + step)
        }
        labels <- breaks[-length(breaks)] + (step / 2)
        levels <- cut(
          raw_values,
          breaks = breaks,
          labels = labels,
          include.lowest = TRUE
        )
      } else {
        levels <- factor(raw_values)
      }

      coeffs <- preds %>%
        mutate(row_id = row_number()) %>%
        group_by(term = !!levels) %>%
        summarise(
          coef = mean(!!fit_colname),
          count = n(),
          se = {
            # Calculate the average design matrix row for this bin
            X_bin_avg <- colMeans(Xi_centered[row_id, , drop = FALSE])
            # Project the covariance matrix: sqrt(avg_row %*% V %*% avg_row_T)
            sqrt(t(X_bin_avg) %*% Vi %*% X_bin_avg)
          },
          .groups = "drop"
        ) %>%
        mutate(
          lower = coef - (1.96 * se),
          upper = coef + (1.96 * se)
        )

      # Reorder levels according to coefficients for factors
      if (any(sapply(preds[term_stripped], is.factor)) && !(grepl('month', term_label))) {
        coeffs <- coeffs %>%
          arrange(coef)

        coeffs$term <- factor(coeffs$term, levels = coeffs$term, ordered = TRUE)

        levels <- factor(levels, levels = coeffs$term, ordered = TRUE)

        # ---------- special case: interaction including month --------------
        # Sort dataframe by 'other' then by month if there is a month interaction term
      } else if (length(term_stripped) == 2 && grepl('month', term_label)) {
        month_idx <- grep('month', term_stripped)
        coeffs <- coeffs %>%
          tidyr::separate(term, into = c("col1", "col2"), sep = ":", remove = FALSE) %>%
          mutate(
            # Assign the month and 'other' based on the original column order
            month_str = if (month_idx == 1) col1 else col2,
            other_str = if (month_idx == 1) col2 else col1,
            # Assign levels to month factor
            month_str = factor(month_str, levels = levels(preds$month), ordered = TRUE)            
            ) %>%
          arrange(other_str, month_str) %>%
          select(-col1, -col2, -month_str, -other_str)
        coeffs$term <- factor(coeffs$term, levels = unique(coeffs$term))
        }
      
      # ------ Organise comparison data for this term if exist ------

      if (
        compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
      ) {
        # Extract components from the comparison list

        if(length(term_stripped)==2) {
        comp_raw_values <- paste(compare_preds_df[[ term_stripped[1] ]], 
                   compare_preds_df[[ term_stripped[2] ]], 
                   sep = ":")
          } else {
          comp_raw_values <- compare_preds_df[[term_stripped]]
          }
                
        comp_preds <- compare_preds_list$preds
        comp_V <- compare_preds_list$V
        comp_X_centered <- compare_preds_list$X_centered
        comp_assign <- compare_preds_list$assign
        comp_terms <- compare_preds_list$terms

        # Numeric terms are cut into factors
        if (
          is.numeric(comp_raw_values) &&
            (!all(comp_raw_values %% 1 == 0) |
              length(unique(comp_raw_values)) > 10)
        ) {
          comp_levels <- cut(
            comp_raw_values,
            breaks = breaks,
            labels = labels,
            include.lowest = TRUE
          )
        } else {
          comp_levels <- factor(comp_raw_values)
        }

        # Identify columns in the comparison model matrix
        comp_match_idx <- which(comp_terms == term_label)
        comp_cols <- which(comp_assign == comp_match_idx)
        # Extract the sub model matrix for this term
        comp_Xi_centered <- comp_X_centered[, comp_cols, drop = FALSE]
        comp_Vi <- comp_V[comp_cols, comp_cols, drop = FALSE]

        comp_coeffs <- compare_preds_df %>%
          mutate(row_id = row_number()) %>%
          group_by(term = !!comp_levels) %>%
          summarise(
            coef = mean(!!fit_colname),
            count = n(),
            se = {
              X_bin_avg <- colMeans(comp_Xi_centered[row_id, , drop = FALSE])
              as.numeric(sqrt(t(X_bin_avg) %*% comp_Vi %*% X_bin_avg))
            },
            .groups = "drop"
          ) %>%
          mutate(
            term = factor(term, levels = levels(coeffs$term), ordered = TRUE),
            lower = coef - (1.96 * se),
            upper = coef + (1.96 * se)
          )
        
        # ------------- Populate Traffic Light dataframe ------------

        # Generate coeffsDiffer indicator
        coeffs_merged <- merge(
          coeffs,
          comp_coeffs,
          by = "term",
          suffixes = c("_new", "_old")
        )

        if (is.factor(raw_values)) {
          coeffsDiffer <- with(
            coeffs_merged,
            abs((exp(coef_new) - exp(coef_old)) / exp(coef_old) * 100)
          )
        } else {
          # weighted averages
          w_new <- with(
            coeffs_merged,
            sum(exp(coef_new) * count_new) / sum(count_new)
          )
          w_old <- with(
            coeffs_merged,
            sum(exp(coef_old) * count_old) / sum(count_old)
          )
          # per cent difference in weighted averages
          coeffsDiffer <- abs((w_new - w_old) / w_old) * 100
        }

        indicators <- data.frame(
          coeffsDiffer = case_when(
            any(coeffsDiffer > 20) ~ '!',
            any(coeffsDiffer > 10) ~ '?',
            TRUE ~ '✔'
          )
        )
      }

      # (1) Coefficients  Plot
      #_____________________________________________________________________________

      # Colours when comparison, black otherwise?
      col1 <- if (compareOn) 'dodgerblue1' else 'black'

      # Create plot theme depending on the term:

      if ((length(levels(levels)) > 12)) {
        # Vertical labels and tighter margins
        dynamic_theme <- theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.x.top = element_text(vjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 15))
        )
      } else {
        # Horizontal labels
        dynamic_theme <- theme(
          axis.text.x = element_text(angle = 0),
          axis.title.x = element_text(margin = margin(t = 10))
        )
      }

      # Create theme elements used for each of the 3 plots
      common_theme <- theme(
        axis.title.x = element_text(vjust = 0),
        panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
        axis.line.x = element_line(colour = "black", linewidth = 0.6),
        axis.line.y = element_line(colour = "black", linewidth = 0.6)
      )

      # No x-labels for vessel_key
      x_labels <- if (grepl('vessel_key', term_label, ignore.case = TRUE)) {
        NULL
        # }  else if(length(levels(levels)) > 15) {
        #   # Downsample labels for continuous vars otherwise they overlap, the logic is to keep every n-th label where n = original no of levels/10
        #   function(x) {
        #     x[seq_along(x) %% max(1, round(length(x) / 10)) != 1] <- ""
        #     return(x)
        #   }
      } else {
        waiver()
      }

      p1 <- ggplot(coeffs, aes(x = term, y = exp(coef))) +
        geom_point(shape = 2, size = 2, color = col1) +
        geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), width = 0, color = col1) +
        # Bottom Cap
        geom_segment(aes(
          x = as.numeric(term) - 0.05,
          xend = as.numeric(term) + 0.05,
          y = exp(lower),
          yend = exp(lower)
        ), color = col1) +
        # Top Cap
        geom_segment(aes(
          x = as.numeric(term) - 0.05,
          xend = as.numeric(term) + 0.05,
          y = exp(upper),
          yend = exp(upper)
        ), color = col1) +
        geom_hline(yintercept = 1, linetype = "dashed") +

        # Conditional block to display coeffs for a model to be compared with (e.g. last year)
        {
          if (
            compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
          ) {
            list(
              geom_point(
                data = comp_coeffs,
                aes(x = term, y = exp(coef)),
                shape = 17,
                size = 2,
                color = '#E41A1CCC'
              ),
              geom_errorbar(
                data = comp_coeffs,
                aes(ymin = exp(lower), ymax = exp(upper)),
                color = '#E41A1CCC',
                width = 0
              )
            )
          }
        } +

        # Trying to place enough labels on the y-scale
        # Solution 1:
        scale_y_log10(breaks = scales::breaks_extended(n = 9)) +

        # Solution 2:
        # {
        #   if (max(exp(coeffs$coef)) - min(exp(coeffs$coef)) > 1) {
        #     scale_y_log10(breaks = scales::breaks_log(n = 8, base = 1.2),
        #                   labels = scales::label_number(accuracy = 0.01))
        #   } else {
        #     scale_y_log10(n.breaks = 8)
        #   }
        # } +
        scale_x_discrete(drop = FALSE, labels = x_labels, position = 'top') +
        theme_bw() +
        background_grid(major = "x", minor = "none") +
        labs(y = 'Effect', x = NULL) +
        dynamic_theme +
        common_theme

      # (2) Distribution Plot
      #_____________________________________________________________________________

      distrs <- preds %>%
        group_by(focus = .[[year]], term = levels) %>%
        summarise(count = n(), .groups = "drop_last") %>%
        mutate(
          total = sum(count),
          prop = count / total
        ) %>%
        ungroup()

      # ---------- special case: interaction including month --------------

      if (length(term_stripped) == 2 && grepl('month', term_label)) {
        
        distrs$term <- factor(distrs$term, levels = unique(coeffs$term))
        
        }

      if (
        compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
      ) {
        comp_distrs <- comp_preds %>%
          group_by(focus = .[[year]], term = comp_levels) %>%
          summarise(count = n(), .groups = "drop_last") %>%

          mutate(
            total = sum(count),
            prop = count / total
          ) %>%

          ungroup()
      }

      p2 <- ggplot(distrs, aes(x = term, y = focus, size = sqrt(prop) * 20)) +
        geom_point(pch = 1, col = col1) +

        # Conditional block to display dists for a model to be compared with (e.g. last year)
        {
          if (
            compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
          ) {
            geom_point(
              data = comp_distrs,
              pch = 1,
              color = '#E41A1CCC',
              alpha = 0.6
            )
          }
        } +

        scale_size_identity() +
        scale_x_discrete(drop = FALSE, labels = x_labels) +
        theme_bw() +
        background_grid(major = "xy", minor = "none") +
        labs(x = paste(term_stripped, collapse = ":"), y = year) +
        dynamic_theme +
        common_theme

      # (3) Influence
      #______________________________________________________________________________

      # This is formula 3b from Bentley et al. 2012
      infl <- preds %>%
        group_by(level = !!sym(year)) %>%
        summarise(Influence = mean(!!sym(fit_colname)))

      # Mean absolute magnitude of the term's impact, converted to annual % change
      overall <- exp(mean(abs(infl$Influence))) - 1

      # Annual growth % of this term's contribution
      n <- nrow(infl)
      trend <- exp(cov(1:n, infl$Influence) / var(1:n)) - 1

      infl <- infl %>%
        mutate(
          Influence = exp(Influence),
          Term = term_label,
          overall = overall,
          trend = trend
        ) %>%
        relocate(Term, .after = level)

      if (
        compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
      ) {
        comp_infl <- compare_preds_df %>%
          group_by(level = !!sym(year)) %>%
          summarise(Influence = exp(mean(!!sym(fit_colname))))
        comp_x_min <- min(comp_infl$Influence)
        comp_x_max <- max(comp_infl$Influence)

        # Generate sign flip indicator:
        influAdded <- infl$Influence[(nrow(comp_infl) + 1):nrow(infl)]
        signLast <- sign(1 - last(comp_infl$Influence))
        signNew <- sign(1 - influAdded)
        indicators$signFlips <- any(signLast != signNew)

        # Generate outofRange indicator
        dist_above <- pmax(infl$Influence - max(comp_infl$Influence), 0)
        dist_below <- pmax(min(comp_infl$Influence) - infl$Influence, 0)
        dist_toRange <- max(dist_above + dist_below)
        indicators$outofRange <- case_when(
          dist_toRange <= 0.05 ~ '✔',
          dist_toRange > 0.1 ~ '!',
          dist_toRange <= 0.1 ~ '?'
        )
      }

      p3 <- ggplot(infl, aes(x = Influence, y = level)) +

        # Conditional block to display coeffs for a model to be compared with (e.g. last year)
        # In the background: The two rectangles (if compareOn is TRUE)
        {
          if (
            compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
          ) {
            list(
              # Outer shaded rectangle (+/- 0.05 padding)
              annotate(
                "rect",
                xmin = comp_x_min - 0.05,
                xmax = comp_x_max + 0.05,
                ymin = -Inf,
                ymax = Inf,
                fill = "grey90",
                alpha = 0.5
              ),
              # Inner shaded rectangle (exact range)
              annotate(
                "rect",
                xmin = comp_x_min,
                xmax = comp_x_max,
                ymin = -Inf,
                ymax = Inf,
                fill = "grey80",
                alpha = 0.5
              )
            )
          }
        } +

        # In the middle layer: Black dots and lines for current (or only) model
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_line(group = 1, col = col1) +
        geom_point(size = 3, pch = 16, col = col1) +

        # In the foreground: Brown dots and lines for model to compare with
        {
          if (
            compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
          ) {
            list(
              geom_line(
                data = comp_infl,
                group = 1,
                linetype = 'dashed',
                color = '#E41A1CCC',
                alpha = 0.6
              ),
              geom_point(
                data = comp_infl,
                shape = 16,
                size = 3,
                color = '#E41A1CCC',
                alpha = 0.6
              )
            )
          }
        } +

        scale_x_continuous(
          limits = function(x) {
            c(min(pretty(x), 0.8), max(pretty(x), 1.2))
          },
          labels = scales::label_number(accuracy = 0.01)
        ) +
        scale_y_discrete(position = "right") +
        labs(x = "Influence", y = NULL) +
        theme_bw() +
        background_grid(major = "y", minor = "none") +
        common_theme

      # (4) Blank or legend for top right corner
      #______________________________________________________________________________
      pv <- ggplot() +
        geom_point(aes(x = 0, y = 0), alpha = 0) +

        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_line(colour = "black", linewidth = 0.6),
          axis.line.y = element_line(colour = "black", linewidth = 0.6)
        )

      # Legend for comparison plot

      if (
        compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
      ) {
        labels = c(
          "Current",
          paste(
            "Previous update to",
            max(as.numeric(as.character(comp_infl$level)))
          ),
          "Range of previous update",
          "Range of \n previous update +5%"
        )

        pv <- ggplot() +
          # Row 1: Current
          annotate(
            "segment",
            x = 5,
            xend = 15,
            y = 50,
            yend = 50,
            color = "dodgerblue1",
            linetype = "solid"
          ) +
          annotate(
            "point",
            x = 10,
            y = 50,
            shape = 16,
            size = 2,
            color = "dodgerblue1"
          ) +
          annotate(
            "point",
            x = 3,
            y = 50,
            shape = 2,
            size = 2,
            color = "dodgerblue1"
          ) +

          # Row 2: Previous
          annotate(
            "segment",
            x = 5,
            xend = 15,
            y = 40,
            yend = 40,
            color = "#E41A1CCC",
            linetype = "dotted"
          ) +
          annotate(
            "point",
            x = 10,
            y = 40,
            shape = 16,
            size = 2,
            color = "#E41A1CCC"
          ) +
          annotate(
            "point",
            x = 3,
            y = 40,
            shape = 17,
            size = 2,
            color = "#E41A1CCC"
          ) +

          # Row 3 & 4: Ranges
          annotate(
            "point",
            x = 10,
            y = 30,
            shape = 15,
            size = 6,
            color = "grey80",
            alpha = 0.5
          ) +
          annotate(
            "point",
            x = 10,
            y = 20,
            shape = 15,
            size = 6,
            color = "grey90",
            alpha = 0.5
          ) +

          # Text Labels
          annotate(
            "text",
            x = 16,
            y = seq(from = 50, to = 20, by = -10),
            label = labels,
            hjust = 0,
            size = 4
          ) +

          coord_cartesian(xlim = c(0, 60), ylim = c(0, 60), expand = FALSE) +

          theme_minimal() +
          theme(
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line.x = element_line(colour = "black", linewidth = 0.6),
            axis.line.y = element_line(colour = "black", linewidth = 0.6)
          )
      }

      combined_plot <- p1 +
        pv +
        p2 +
        p3 +
        plot_layout(nrow = 2, ncol = 2, heights = c(1, 2), widths = c(2, 1)) &
        theme(
          # axis.title = element_text(size = 14),
          #  axis.text = element_text(size = 14),
          plot.margin = margin(0, 0, 0, 0)
        )

      # print(combined_plot)
      attr(combined_plot, "indicators") <- indicators
      combined_plot@meta$indicators <- indicators
      return(combined_plot)
    }),
    terms_labels
  )

  return(result)
}


# Draft:
#   cdi_plot_with_indicators <- function(preds_list, compare_preds_list = NULL) {
#   # ------------------- Prepare the data ----------------------
#   # comparison switch
#   compareOn <- !is.null(compare_preds_list)

#   if (compareOn) {
#     compare_preds_df <- compare_preds_list$preds
#   }

#   indicators <- NULL # initialise indicators dataframe populateted if compareOn

#   # Extract components from the preds list

#   preds <- preds_list$preds
#   V <- preds_list$V
#   X_centered <- preds_list$X_centered
#   assigns <- preds_list$assign
#   all_model_terms <- preds_list$terms
#   year <- preds_list$year

#   # exclude year from terms
#   terms_labels <- setdiff(all_model_terms, year)

#   # ---------- Loop through terms to generate plots ------------

#   result <- setNames(
#     lapply(terms_labels, function(term_label) {
#       # ---------- Housekeeping  ------------------
#       # Define col names for fitted terms and standard errors
#       fit_colname <- sym(paste0('fit.', term_label))
#       se_colname <- sym(paste0('se.fit.', term_label))

#       # Identify model matrix columns for this specific term
#       match_idx <- which(all_model_terms == term_label)
#       cols <- which(assigns == match_idx)

#       # Extract the sub model matrix for this term
#       Xi_centered <- X_centered[, cols, drop = FALSE]
#       Vi <- V[cols, cols, drop = FALSE]

#       # Define levels of term on which coefficient and distribution plots will be based
#       # This is done by search for each column name in the data as a whole word in the
#       # each term. This allows for matching of terms like 'poly(log(var),3)' with 'var'

#       # ---------- Assess term type ----------------------------

#       term_stripped <- all.vars(as.formula(paste("~", term_label)))
#       is_interaction <- length(term_stripped) == 2
#       is_numeric_interaction <- is_interaction &&
#         !all(sapply(preds[term_stripped], is.factor))

#       # ---------- Organise data for this term ------------------

#       # Helper function to bin numeric vectors
#       bin_numeric <- function(x, n_breaks) {
#         if (is.numeric(x) && (!all(x %% 1 == 0) | length(unique(x)) > 10)) {
#           breaks <- pretty(x, n_breaks)
#           step <- breaks[2] - breaks[1]
#           breaks <- breaks[breaks <= max(x)]
#           if (max(breaks) < max(x)) {
#             breaks <- c(breaks, max(breaks) + step)
#           }
#           labels <- breaks[-length(breaks)] + (step / 2)
#           return(cut(
#             x,
#             breaks = breaks,
#             labels = labels,
#             include.lowest = TRUE
#           ))
#         }
#         return(factor(x))
#       }

#       if (is_numeric_interaction) {
#         # Interaction Case when at least one interaction term is numeric
#         first_numeric_term <- Find(
#           function(x) !is.factor(preds[[x]]),
#           term_stripped
#         )
#         second_term <- Find(function(x) x != first_numeric_term, term_stripped)
#         raw_values_1 <- preds[[first_numeric_term]]
#         raw_values_2 <- preds[[second_term]]

#         # A) Slice first numeric term into Quantiles (Q1, Median/IQR, Q3) to act as groups
#         quant_breaks <- quantile(
#           raw_values_1,
#           probs = c(0, 0.25, 0.75, 1),
#           na.rm = TRUE
#         )

#         levels_1 <- cut(
#           raw_values_1,
#           breaks = quant_breaks,
#           include.lowest = TRUE
#         )
#         # Replace default labels with medians
#         bin_medians <- tapply(raw_values_1, levels_1, median)
#         levels(levels_1) <- bin_medians

#         # B) Handle second term normally
#         # Numeric terms are cut into factors
#       if (
#         is.numeric(raw_values_2) &&
#           (!all(raw_values_2 %% 1 == 0) | length(unique(raw_values_2)) > 10)
#       ) {
#         breaks <- pretty(raw_values_2, 15)
#         step <- breaks[2] - breaks[1]
#         breaks <- breaks[breaks <= max(raw_values_2)]

#         if (max(breaks) < max(raw_values_2)) {
#           breaks <- c(breaks, max(breaks) + step)
#         }
#         labels <- breaks[-length(breaks)] + (step / 2)
#         levels_2 <- cut(
#           raw_values_2,
#           breaks = breaks,
#           labels = labels,
#           include.lowest = TRUE
#         )
#       } else {
#         levels_2 <- factor(raw_values_2)
#       }
        
#         # Paste levels of both predicoturs for using in the plots
#         levels <- factor(paste(levels_1, levels_2, sep = ":"))
#       } else {
#         if (is_interaction) {
#           # Interaction Case when both interaction terms are factors
#           raw_values <- paste(
#             preds[[term_stripped[1]]],
#             preds[[term_stripped[2]]],
#             sep = ":"
#           )
#         } else {
#           # Not an interaction case
#           raw_values <- preds[[term_stripped]]
#         }
#         # Numeric terms are cut into factors
#       if (
#         is.numeric(raw_values) &&
#           (!all(raw_values %% 1 == 0) | length(unique(raw_values)) > 10)
#       ) {
#         breaks <- pretty(raw_values, 30)
#         step <- breaks[2] - breaks[1]
#         breaks <- breaks[breaks <= max(raw_values)]

#         if (max(breaks) < max(raw_values)) {
#           breaks <- c(breaks, max(breaks) + step)
#         }
#         labels <- breaks[-length(breaks)] + (step / 2)
#         levels <- cut(
#           raw_values,
#           breaks = breaks,
#           labels = labels,
#           include.lowest = TRUE
#         )
#       } else {
#         levels <- factor(raw_values)
#       }
#       }

#       coeffs <- preds %>%
#         mutate(row_id = row_number()) %>%
#         group_by(term = !!levels) %>%
#         summarise(
#           coef = mean(!!fit_colname),
#           count = n(),
#           se = {
#             # Calculate the average design matrix row for this bin
#             X_bin_avg <- colMeans(Xi_centered[row_id, , drop = FALSE])
#             # Project the covariance matrix: sqrt(avg_row %*% V %*% avg_row_T)
#             sqrt(t(X_bin_avg) %*% Vi %*% X_bin_avg)
#           },
#           .groups = "drop"
#         ) %>%
#         mutate(
#           lower = coef - (1.96 * se),
#           upper = coef + (1.96 * se)
#         )

#       # ----- Reorder levels according to coefficients ------------
#       if (sum(sapply(preds[term_stripped], is.factor)) == length(term_stripped) && !(grepl('month', term_label))
#       ) {

#        # --- A). Factor, or factor:factor, excluding 'month' ---
#         # Sorts by the effect size 
#         coeffs <- coeffs %>%
#           arrange(coef)
                
#       } else if (is_interaction && grepl('month', term_label)) {
#         # --- B). interaction including month ---
#         # Sort dataframe by 'other' then by month if there is a month interaction term
#         month_idx <- grep('month', term_stripped)
#         coeffs <- coeffs %>%
#           tidyr::separate(
#             term,
#             into = c("col1", "col2"),
#             sep = ":",
#             remove = FALSE
#           ) %>%
#           mutate(
#             # Assign the month and 'other' based on the original column order
#             month_str = if (month_idx == 1) col1 else col2,
#             other_str = if (month_idx == 1) col2 else col1,
#             # Assign levels to month factor
#             month_str = factor(
#               month_str,
#               levels = levels(preds$month),
#               ordered = TRUE
#             )
#           ) %>%
#           arrange(other_str, month_str) %>%
#           select(-col1, -col2, -month_str, -other_str)
        
#       } else if (is_numeric_interaction) {
#         # --- C). numeric:numeric interaction ---
#         coeffs <- coeffs %>%
#           tidyr::separate(term, into = c("num1", "num2"), sep = ":", remove = FALSE) %>%
#           mutate(
#             num1 = suppressWarnings(as.numeric(num1)),
#             num2 = suppressWarnings(as.numeric(num2))
#           ) %>%
#           arrange(num1, num2) %>%
#           select(-num1, -num2)
          
        
#       } else if (is_interaction ) {
#         # --- D). factor:numeric interaction ---
#         # Find which position holds the factor and which holds the numeric
#         factor_idx <- which(sapply(preds[term_stripped], is.factor))
#         num_idx <- which(!sapply(preds[term_stripped], is.factor))

#         coeffs <- coeffs %>%
#           tidyr::separate(term, into = c("val1", "val2"), sep = ":", remove = FALSE) %>%
#           mutate(
#             fact_str = if (factor_idx == 1) val1 else val2,
#             num_str  = if (num_idx == 1) val1 else val2,
            
#             # Apply original factor levels so it sorts correctly
#             fact_str = factor(fact_str, levels = levels(preds[[ term_stripped[factor_idx] ]])),
#             num_str  = suppressWarnings(as.numeric(num_str))
#           ) %>%
#           # Arrange by numeric values, then by factor
#           arrange(num_val, fact_str) %>%
#           select(-val1, -val2, -fact_str, -num_str)
          
        
#       }
#       coeffs$term <- factor(coeffs$term, levels = unique(coeffs$term), ordered = TRUE)
#       levels <- factor(levels, levels = coeffs$term, ordered = TRUE)


#       # ------ Organise comparison data for this term if exist ------

#       if (
#         compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
#       ) {
#         # Extract components from the comparison list

#         if (is_numeric_interaction) {
#            raw_values_1 <- compare_preds_df[[first_numeric_term]]
#         raw_values_2 <- compare_preds_df[[second_term]]
#           # A) Slice first numeric term into Quantiles (Q1, Median/IQR, Q3) to act as groups.
#           # # Important: ise the same qualtile breaks as the first model

#         levels_1 <- cut(
#           raw_values_1,
#           breaks = quant_breaks,
#           labels = bin_medians,
#           include.lowest = TRUE
#         )
        
#         # B) Handle second term normally

#         # Numeric terms are cut into factors
#         if (
#           is.numeric(raw_values_2) &&
#             (!all(raw_values_2 %% 1 == 0) |
#               length(unique(raw_values_2)) > 10)
#         ) {
#           levels_2 <- cut(
#             raw_values_2,
#             breaks = breaks,
#             labels = labels,
#             include.lowest = TRUE
#           )
#         } else {
#           levels_2 <- factor(raw_values_2)
#         }

#         # Paste levels of both predicoturs for using in the plots
#         comp_levels <- factor(paste(levels_1, levels_2, sep = ":"))
#         } else {
#           if (is_interaction) {
#           # Interaction Case when both interaction terms are factors
#           comp_raw_values <- paste(
#             compare_preds_df[[term_stripped[1]]],
#             compare_preds_df[[term_stripped[2]]],
#             sep = ":"
#           )
#         } else {
#           # Not an interaction case
#           comp_raw_values <- compare_preds_df[[term_stripped]]
#         }
#          # Numeric terms are cut into factors
#         if (
#           is.numeric(comp_raw_values) &&
#             (!all(comp_raw_values %% 1 == 0) |
#               length(unique(comp_raw_values)) > 10)
#         ) {
#           comp_levels <- cut(
#             comp_raw_values,
#             breaks = breaks,
#             labels = labels,
#             include.lowest = TRUE
#           )
#         } else {
#           comp_levels <- factor(comp_raw_values)
#         }
#         }

               
#         comp_preds <- compare_preds_list$preds
#         comp_V <- compare_preds_list$V
#         comp_X_centered <- compare_preds_list$X_centered
#         comp_assign <- compare_preds_list$assign
#         comp_terms <- compare_preds_list$terms

#         # Identify columns in the comparison model matrix
#         comp_match_idx <- which(comp_terms == term_label)
#         comp_cols <- which(comp_assign == comp_match_idx)
#         # Extract the sub model matrix for this term
#         comp_Xi_centered <- comp_X_centered[, comp_cols, drop = FALSE]
#         comp_Vi <- comp_V[comp_cols, comp_cols, drop = FALSE]

#         comp_coeffs <- compare_preds_df %>%
#           mutate(row_id = row_number()) %>%
#           group_by(term = !!comp_levels) %>%
#           summarise(
#             coef = mean(!!fit_colname),
#             count = n(),
#             se = {
#               X_bin_avg <- colMeans(comp_Xi_centered[row_id, , drop = FALSE])
#               as.numeric(sqrt(t(X_bin_avg) %*% comp_Vi %*% X_bin_avg))
#             },
#             .groups = "drop"
#           ) %>%
#           mutate(
#             term = factor(term, levels = levels(coeffs$term), ordered = TRUE),
#             lower = coef - (1.96 * se),
#             upper = coef + (1.96 * se)
#           )

#         # ------------- Populate Traffic Light dataframe ------------

#         # Generate coeffsDiffer indicator
#         coeffs_merged <- merge(
#           coeffs,
#           comp_coeffs,
#           by = "term",
#           suffixes = c("_new", "_old")
#         )

#         if (any(sapply(preds[term_stripped], is.factor))) {
#           coeffsDiffer <- with(
#             coeffs_merged,
#             abs((exp(coef_new) - exp(coef_old)) / exp(coef_old) * 100)
#           )
#         } else {
#           # weighted averages
#           w_new <- with(
#             coeffs_merged,
#             sum(exp(coef_new) * count_new) / sum(count_new)
#           )
#           w_old <- with(
#             coeffs_merged,
#             sum(exp(coef_old) * count_old) / sum(count_old)
#           )
#           # per cent difference in weighted averages
#           coeffsDiffer <- abs((w_new - w_old) / w_old) * 100
#         }

#         indicators <- data.frame(
#           coeffsDiffer = case_when(
#             any(coeffsDiffer > 20) ~ '!',
#             any(coeffsDiffer > 10) ~ '?',
#             TRUE ~ '✔'
#           )
#         )
#         }

#       # (1) Coefficients  Plot
#       #_____________________________________________________________________________

#       # Colours when comparison, black otherwise?
#       col1 <- if (compareOn) 'dodgerblue1' else 'black'

#       # Create plot theme depending on the term:

#       if ((length(levels(levels)) > 12)) {
#         # Vertical labels and tighter margins
#         dynamic_theme <- theme(
#           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#           axis.text.x.top = element_text(vjust = 0.5),
#           axis.title.x = element_text(margin = margin(t = 15))
#         )
#       } else {
#         # Horizontal labels
#         dynamic_theme <- theme(
#           axis.text.x = element_text(angle = 0),
#           axis.title.x = element_text(margin = margin(t = 10))
#         )
#       }

#       # Create theme elements used for each of the 3 plots
#       common_theme <- theme(
#         axis.title.x = element_text(vjust = 0),
#         panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
#         axis.line.x = element_line(colour = "black", linewidth = 0.6),
#         axis.line.y = element_line(colour = "black", linewidth = 0.6)
#       )

#       # No x-labels for vessel_key
#       x_labels <- if (grepl('vessel_key', term_label, ignore.case = TRUE)) {
#         NULL
#         # }  else if(length(levels(levels)) > 15) {
#         #   # Downsample labels for continuous vars otherwise they overlap, the logic is to keep every n-th label where n = original no of levels/10
#         #   function(x) {
#         #     x[seq_along(x) %% max(1, round(length(x) / 10)) != 1] <- ""
#         #     return(x)
#         #   }
#       } else {
#         waiver()
#       }

#       p1 <- ggplot(coeffs, aes(x = term, y = exp(coef))) +
#         geom_point(shape = 2, size = 2, color = col1) +
#         geom_errorbar(
#           aes(ymin = exp(lower), ymax = exp(upper)),
#           width = 0,
#           color = col1
#         ) +
#         # Bottom Cap
#         geom_segment(
#           aes(
#             x = as.numeric(term) - 0.05,
#             xend = as.numeric(term) + 0.05,
#             y = exp(lower),
#             yend = exp(lower)
#           ),
#           color = col1
#         ) +
#         # Top Cap
#         geom_segment(
#           aes(
#             x = as.numeric(term) - 0.05,
#             xend = as.numeric(term) + 0.05,
#             y = exp(upper),
#             yend = exp(upper)
#           ),
#           color = col1
#         ) +
#         geom_hline(yintercept = 1, linetype = "dashed") +

#         # Conditional block to display coeffs for a model to be compared with (e.g. last year)
#         {
#           if (
#             compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
#           ) {
#             list(
#               geom_point(
#                 data = comp_coeffs,
#                 aes(x = term, y = exp(coef)),
#                 shape = 17,
#                 size = 2,
#                 color = '#E41A1CCC'
#               ),
#               geom_errorbar(
#                 data = comp_coeffs,
#                 aes(ymin = exp(lower), ymax = exp(upper)),
#                 color = '#E41A1CCC',
#                 width = 0
#               )
#             )
#           }
#         } +

#         # Trying to place enough labels on the y-scale
#         # Solution 1:
#         scale_y_log10(breaks = scales::breaks_extended(n = 9)) +

#         # Solution 2:
#         # {
#         #   if (max(exp(coeffs$coef)) - min(exp(coeffs$coef)) > 1) {
#         #     scale_y_log10(breaks = scales::breaks_log(n = 8, base = 1.2),
#         #                   labels = scales::label_number(accuracy = 0.01))
#         #   } else {
#         #     scale_y_log10(n.breaks = 8)
#         #   }
#         # } +
#         scale_x_discrete(drop = FALSE, labels = x_labels, position = 'top') +
#         theme_bw() +
#         background_grid(major = "x", minor = "none") +
#         labs(y = 'Effect', x = NULL) +
#         dynamic_theme +
#         common_theme

#       # (2) Distribution Plot
#       #_____________________________________________________________________________

#       distrs <- preds %>%
#         group_by(focus = .[[year]], term = levels) %>%
#         summarise(count = n(), .groups = "drop_last") %>%
#         mutate(
#           total = sum(count),
#           prop = count / total
#         ) %>%
#         ungroup()

#       # ---------- special case: interaction including month --------------

#       if (length(term_stripped) == 2 && grepl('month', term_label)) {
#         distrs$term <- factor(distrs$term, levels = unique(coeffs$term))
#       }

#       if (
#         compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
#       ) {
#         comp_distrs <- comp_preds %>%
#           group_by(focus = .[[year]], term = comp_levels) %>%
#           summarise(count = n(), .groups = "drop_last") %>%

#           mutate(
#             total = sum(count),
#             prop = count / total
#           ) %>%

#           ungroup()
#       }

#       p2 <- ggplot(distrs, aes(x = term, y = focus, size = sqrt(prop) * 20)) +
#         geom_point(pch = 1, col = col1) +

#         # Conditional block to display dists for a model to be compared with (e.g. last year)
#         {
#           if (
#             compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
#           ) {
#             geom_point(
#               data = comp_distrs,
#               pch = 1,
#               color = '#E41A1CCC',
#               alpha = 0.6
#             )
#           }
#         } +

#         scale_size_identity() +
#         scale_x_discrete(drop = FALSE, labels = x_labels) +
#         theme_bw() +
#         background_grid(major = "xy", minor = "none") +
#         labs(x = paste(term_stripped, collapse = ":"), y = year) +
#         dynamic_theme +
#         common_theme

#       # (3) Influence
#       #______________________________________________________________________________

#       # This is formula 3b from Bentley et al. 2012
#       infl <- preds %>%
#         group_by(level = !!sym(year)) %>%
#         summarise(Influence = mean(!!sym(fit_colname)))

#       # Mean absolute magnitude of the term's impact, converted to annual % change
#       overall <- exp(mean(abs(infl$Influence))) - 1

#       # Annual growth % of this term's contribution
#       n <- nrow(infl)
#       trend <- exp(cov(1:n, infl$Influence) / var(1:n)) - 1

#       infl <- infl %>%
#         mutate(
#           Influence = exp(Influence),
#           Term = term_label,
#           overall = overall,
#           trend = trend
#         ) %>%
#         relocate(Term, .after = level)

#       if (
#         compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
#       ) {
#         comp_infl <- compare_preds_df %>%
#           group_by(level = !!sym(year)) %>%
#           summarise(Influence = exp(mean(!!sym(fit_colname))))
#         comp_x_min <- min(comp_infl$Influence)
#         comp_x_max <- max(comp_infl$Influence)

#         # Generate sign flip indicator:
#         influAdded <- infl$Influence[(nrow(comp_infl) + 1):nrow(infl)]
#         signLast <- sign(1 - last(comp_infl$Influence))
#         signNew <- sign(1 - influAdded)
#         indicators$signFlips <- any(signLast != signNew)

#         # Generate outofRange indicator
#         dist_above <- pmax(infl$Influence - max(comp_infl$Influence), 0)
#         dist_below <- pmax(min(comp_infl$Influence) - infl$Influence, 0)
#         dist_toRange <- max(dist_above + dist_below)
#         indicators$outofRange <- case_when(
#           dist_toRange <= 0.05 ~ '✔',
#           dist_toRange > 0.1 ~ '!',
#           dist_toRange <= 0.1 ~ '?'
#         )
#       }

#       p3 <- ggplot(infl, aes(x = Influence, y = level)) +

#         # Conditional block to display coeffs for a model to be compared with (e.g. last year)
#         # In the background: The two rectangles (if compareOn is TRUE)
#         {
#           if (
#             compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
#           ) {
#             list(
#               # Outer shaded rectangle (+/- 0.05 padding)
#               annotate(
#                 "rect",
#                 xmin = comp_x_min - 0.05,
#                 xmax = comp_x_max + 0.05,
#                 ymin = -Inf,
#                 ymax = Inf,
#                 fill = "grey90",
#                 alpha = 0.5
#               ),
#               # Inner shaded rectangle (exact range)
#               annotate(
#                 "rect",
#                 xmin = comp_x_min,
#                 xmax = comp_x_max,
#                 ymin = -Inf,
#                 ymax = Inf,
#                 fill = "grey80",
#                 alpha = 0.5
#               )
#             )
#           }
#         } +

#         # In the middle layer: Black dots and lines for current (or only) model
#         geom_vline(xintercept = 1, linetype = "dashed") +
#         geom_line(group = 1, col = col1) +
#         geom_point(size = 3, pch = 16, col = col1) +

#         # In the foreground: Brown dots and lines for model to compare with
#         {
#           if (
#             compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
#           ) {
#             list(
#               geom_line(
#                 data = comp_infl,
#                 group = 1,
#                 linetype = 'dashed',
#                 color = '#E41A1CCC',
#                 alpha = 0.6
#               ),
#               geom_point(
#                 data = comp_infl,
#                 shape = 16,
#                 size = 3,
#                 color = '#E41A1CCC',
#                 alpha = 0.6
#               )
#             )
#           }
#         } +

#         scale_x_continuous(
#           limits = function(x) {
#             c(min(pretty(x), 0.8), max(pretty(x), 1.2))
#           },
#           labels = scales::label_number(accuracy = 0.01)
#         ) +
#         scale_y_discrete(position = "right") +
#         labs(x = "Influence", y = NULL) +
#         theme_bw() +
#         background_grid(major = "y", minor = "none") +
#         common_theme

#       # (4) Blank or legend for top right corner
#       #______________________________________________________________________________
#       pv <- ggplot() +
#         geom_point(aes(x = 0, y = 0), alpha = 0) +

#         theme_minimal() +
#         theme(
#           panel.grid = element_blank(),
#           axis.title = element_blank(),
#           axis.text = element_blank(),
#           axis.line.x = element_line(colour = "black", linewidth = 0.6),
#           axis.line.y = element_line(colour = "black", linewidth = 0.6)
#         )

#       # Legend for comparison plot

#       if (
#         compareOn && paste0('fit.', term_label) %in% names(compare_preds_df)
#       ) {
#         labels = c(
#           "Current",
#           paste(
#             "Previous update to",
#             max(as.numeric(as.character(comp_infl$level)))
#           ),
#           "Range of previous update",
#           "Range of \n previous update +5%"
#         )

#         pv <- ggplot() +
#           # Row 1: Current
#           annotate(
#             "segment",
#             x = 5,
#             xend = 15,
#             y = 50,
#             yend = 50,
#             color = "dodgerblue1",
#             linetype = "solid"
#           ) +
#           annotate(
#             "point",
#             x = 10,
#             y = 50,
#             shape = 16,
#             size = 2,
#             color = "dodgerblue1"
#           ) +
#           annotate(
#             "point",
#             x = 3,
#             y = 50,
#             shape = 2,
#             size = 2,
#             color = "dodgerblue1"
#           ) +

#           # Row 2: Previous
#           annotate(
#             "segment",
#             x = 5,
#             xend = 15,
#             y = 40,
#             yend = 40,
#             color = "#E41A1CCC",
#             linetype = "dotted"
#           ) +
#           annotate(
#             "point",
#             x = 10,
#             y = 40,
#             shape = 16,
#             size = 2,
#             color = "#E41A1CCC"
#           ) +
#           annotate(
#             "point",
#             x = 3,
#             y = 40,
#             shape = 17,
#             size = 2,
#             color = "#E41A1CCC"
#           ) +

#           # Row 3 & 4: Ranges
#           annotate(
#             "point",
#             x = 10,
#             y = 30,
#             shape = 15,
#             size = 6,
#             color = "grey80",
#             alpha = 0.5
#           ) +
#           annotate(
#             "point",
#             x = 10,
#             y = 20,
#             shape = 15,
#             size = 6,
#             color = "grey90",
#             alpha = 0.5
#           ) +

#           # Text Labels
#           annotate(
#             "text",
#             x = 16,
#             y = seq(from = 50, to = 20, by = -10),
#             label = labels,
#             hjust = 0,
#             size = 4
#           ) +

#           coord_cartesian(xlim = c(0, 60), ylim = c(0, 60), expand = FALSE) +

#           theme_minimal() +
#           theme(
#             panel.grid = element_blank(),
#             axis.title = element_blank(),
#             axis.text = element_blank(),
#             axis.line.x = element_line(colour = "black", linewidth = 0.6),
#             axis.line.y = element_line(colour = "black", linewidth = 0.6)
#           )
#       }

#       combined_plot <- p1 +
#         pv +
#         p2 +
#         p3 +
#         plot_layout(nrow = 2, ncol = 2, heights = c(1, 2), widths = c(2, 1)) &
#         theme(
#           # axis.title = element_text(size = 14),
#           #  axis.text = element_text(size = 14),
#           plot.margin = margin(0, 0, 0, 0)
#         )

#       # print(combined_plot)
#       attr(combined_plot, "indicators") <- indicators
#       combined_plot@meta$indicators <- indicators
#       return(combined_plot)
#     }),
#     terms_labels
#   )

#   return(result)
# }