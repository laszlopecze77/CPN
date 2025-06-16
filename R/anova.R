#' Analysis of Deviance for CPN Models
#'
#' Computes an analysis of deviance table for objects of class `cpn`, either:
#' 1. sequentially by model terms (Type I ANOVA) when one model is supplied, or
#' 2. by comparing two or more nested models.
#'
#' @param object An object of class `cpn`, typically the result of a call
#'   to [cpn()].
#' @param ... Optional additional objects of class `cpn` for model comparison.
#'
#' @return An object of class `"anova"` inheriting from `"data.frame"`,
#' containing:
#' - If one model is provided: sequential deviance table by added terms.
#' - If multiple models are provided: deviance comparison between models.
#'
#' @details
#' When a single model is supplied, this function refits the model with terms
#' added sequentially and reports the change in residual deviance.
#'
#' When multiple models are supplied, they are assumed to be nested and ordered
#' by increasing complexity (based on residual degrees of freedom). The function
#' then compares their deviances using chi-squared tests.
#'
#' @examples
#' # Simulated data
#' set.seed(123)
#' data <- simulate_cpn_data()
#'
#' # Sequential analysis of deviance
#' fit <- cpn(y ~ x1 + x2, data = data)
#' anova(fit)
#'
#' # Model comparison
#' fit1 <- cpn(y ~ x1, data = data)
#' fit2 <- cpn(y ~ x1 + x2, data = data)
#' anova(fit1, fit2)
#'
#' @export
anova.cpn <- function(object, ...) {
  if (length(list(...)) > 0) {
    # Model comparison
    h_anova_comp(object, ...)
  } else {
    # Sequential analysis of deviance
    h_anova_seq(object)
  }
}

#' Helper function: Analysis of Deviance for Comparison of CPN Models
#'
#' Computes an analysis of deviance table for objects of class `cpn`,
#' by comparing two or more nested models.
#'
#' @param ... `cpn` models to compare.
#'
#' @return An object of class `"anova"` inheriting from `"data.frame"`,
#'   containing sequential deviance table by added terms.
#'
#' @details
#' The multiple supplied models are assumed to be nested and ordered
#' by increasing complexity (based on residual degrees of freedom). The function
#' then compares their deviances using chi-squared tests.
#'
#' @noRd
h_anova_comp <- function(...) {
  models <- list(...)

  # Check all inputs are of class "cpn"
  if (!all(sapply(models, inherits, "cpn"))) {
    stop("All inputs must be of class 'cpn'")
  }

  # Extract and order by residual df (increasing model complexity)
  df_res <- sapply(models, function(m) m$df_residual)
  dev <- sapply(models, function(m) m$residual_deviance)
  n_models <- length(models)

  ord <- order(df_res)
  models <- models[ord]
  df_res <- df_res[ord]
  dev <- dev[ord]

  # Compute differences and p-values
  df_diff <- c(NA, diff(df_res))
  dev_diff <- c(NA, diff(dev))
  p_vals <- c(
    NA,
    stats::pchisq(dev_diff[-1], df_diff[-1], lower.tail = FALSE)
  )

  # Build comparison table
  structure(
    data.frame(
      Model = paste0("Model ", seq_len(n_models)),
      `Resid. Df` = df_res,
      `Resid. Dev` = dev,
      Df = df_diff,
      Deviance = dev_diff,
      `Pr(>Chi)` = p_vals,
      check.names = FALSE
    ),
    class = c("anova.cpn", "data.frame")
  )
}

#' Helper function: Sequential Analysis of Deviance for CPN Models
#'
#' Computes an analysis of deviance table for objects of class `cpn`,
#' sequentially by model terms (Type I ANOVA) when one model is supplied, or
#'
#' @param x a `cpn` model.
#'
#' @return An object of class `"anova"` inheriting from `"data.frame"`,
#'   containing deviance comparison between models.
#'
#' @details
#' This function refits the model with terms added sequentially and reports the
#' change in residual deviance.
#'
#' @noRd
h_anova_seq <- function(x) {

  terms_obj <- stats::terms(x)
  term_labels <- attr(terms_obj, "term.labels")

  if (length(term_labels) == 0) {
    warning("Model contains only an intercept. No terms to test.")
    return(invisible(NULL))
  }

  response <- as.character(attr(terms_obj, "variables"))[2]
  data <- x$data

  # Intercept-only model
  null_formula <- stats::as.formula(paste(response, "~ 1"))
  fit_prev <- cpn(null_formula, data = data)
  dev_prev <- fit_prev$residual_deviance
  df_prev <- fit_prev$df_residual

  rows <- list()
  rows[[1]] <- list(
    Term = "Residuals",
    Df = NA,
    Deviance = NA,
    `Resid. Df` = df_prev,
    `Resid. Dev` = dev_prev,
    `Pr(>Chi)` = NA
  )

  for (i in seq_along(term_labels)) {
    term <- term_labels[i]
    current_formula <- stats::as.formula(
      paste(response, "~", paste(term_labels[1:i], collapse = " + "))
    )
    fit_curr <- cpn(current_formula, data = data)

    dev_curr <- fit_curr$residual_deviance
    df_curr <- fit_curr$df_residual

    df_diff <- df_prev - df_curr
    dev_diff <- dev_prev - dev_curr
    p_val <- if (!is.na(df_diff) && df_diff > 0) {
      stats::pchisq(dev_diff, df = df_diff, lower.tail = FALSE)
    } else {
      NA
    }

    rows[[i + 1]] <- list(
      Term = term,
      Df = df_diff,
      Deviance = dev_diff,
      `Resid. Df` = df_curr,
      `Resid. Dev` = dev_curr,
      `Pr(>Chi)` = p_val
    )

    dev_prev <- dev_curr
    df_prev <- df_curr
  }

  rows_list <- lapply(rows, function(row) {
    as.data.frame(row, check.names = FALSE)
  })

  structure(
    do.call(rbind, rows_list),
    rownames = NULL,
    class = c("anova.cpn", "data.frame")
  )
}

#' @export
print.anova.cpn <- function(x, digits = max(4, getOption("digits") - 2), ...) {
  signif_stars <- if ("Pr(>Chi)" %in% names(x)) {
    stats::symnum(x$`Pr(>Chi)`,
                  corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
  } else {
    NULL
  }

  out <- data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)

  if (!is.null(signif_stars)) {
    out$Signif <- signif_stars
  }

  # Convert NAs to empty strings (only for display)
  out_print <- out
  out_print[] <- lapply(out_print, function(col) {
    if (is.numeric(col)) {
      ifelse(is.na(col), "", formatC(col, digits = digits, format = "g"))
    } else {
      ifelse(is.na(col), "", as.character(col))
    }
  })

  print(out_print, right = TRUE, row.names = FALSE, ...)

  if (!is.null(signif_stars)) {
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }

  invisible(x)
}
