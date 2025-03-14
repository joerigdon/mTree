#' Get averaged matches
#'
#' @param matched_data Filtered matched data set obtained from [`get_filtered_matches()`]
#' @param treatment Character. Name of the treatment column in `data`. Currently
#'    assumes `1` for treatment and `0` for control.
#' @param outcome Character. Name of the outcome column in `data`.
#' @param covariates Character vector. Vector of covariate column names
#'        from `data` to use in match
#' @param match_id Character. The column name in `matched_data` that contains
#'       the matching identifier. Default is "`match_id`".
#' @param weight Character. The column name in `matched_data` that contains
#'       the weight. Default is "`wt`".
#'
#' @return A data frame of averaged matches
#' @export
#'
#' @examples
#' matches <- get_filtered_matches(
#'   example_trial_data,
#'   treatment = "Z",
#'   covariates = c("age", "stage4", "site", "ecog", "prevTrt", "diseaseFree")
#' )
#' get_averaged_matches(
#'   matches,
#'   treatment = "Z",
#'   outcome = "Y",
#'   covariates = c("age", "stage4", "site", "ecog", "prevTrt", "diseaseFree"))
get_averaged_matches <- function(matched_data, treatment, outcome, covariates,
                                 match_id = "match_id", weight = "wt") {
  data <- matched_data
  required_cols <- c(treatment, outcome, match_id, weight)
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    cli::cli_abort(c("x" = "Missing required columns: {missing}"))
  }

  missing_covs <- setdiff(covariates, names(data))
  if (length(missing_covs) > 0) {
    cli::cli_abort(c("x" = "Missing covariates in data: {missing_covs}"))
  }

  unique_matches <- unique(data[[match_id]])
  result <- data.frame(
    match = unique_matches,
    delta = numeric(length(unique_matches)),
    weight = numeric(length(unique_matches))
  )
  names(result)[3] <- weight

  for (cov in covariates) {
    col_class <- class(data[[cov]])
    if ("factor" %in% col_class) {
      result[[cov]] <- factor(rep(NA, length(unique_matches)), levels=levels(data[[cov]]))
    } else {
      result[[cov]] <- vector(mode = typeof(data[[cov]]), length = length(unique_matches))
    }
    class(result[[cov]]) <- col_class
  }
  for (i in seq_along(unique_matches)) {

    match_data <- data[data[[match_id]] == unique_matches[i], ]

    if (sum(match_data[[treatment]] == 1) != 1 || sum(match_data[[treatment]] == 0) != 1) {
      cli::cli_abort(c("x" = "Match {unique_matches[i]} does not have exactly one treatment and one control unit"))
    }

    t_row <- match_data[match_data[[treatment]] == 1, ]
    c_row <- match_data[match_data[[treatment]] == 0, ]

    result$delta[i] <- t_row[[outcome]] - c_row[[outcome]]
    result[[weight]][i] <- t_row[[weight]]

    for (cov in covariates) {
      if (!inherits(data[[cov]],"numeric")) {
        if (t_row[[cov]] != c_row[[cov]]) {
          cli::cli_abort(c("!" = "Categorical variable {cov} has different values for match {unique_matches[i]}"))
        }
        result[[cov]][i] <- t_row[[cov]]
      } else {
        result[[cov]][i] <- (t_row[[cov]] + c_row[[cov]]) / 2
      }
    }
  }

  return(result)
}
