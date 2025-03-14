#' Get filtered matches
#'
#' @param data Data frame.
#' @param treatment Character. Name of the treatment column in `data`. Currently
#'    assumes `1` for treatment and `0` for control.
#' @param covariates Character vector. Vector of covariate column names
#'        from `data` to use in match
#' @param id Character. Name of the id column in `data`
#' @param weight Character. The name that should be given to the variable
#'        containing the weights in the output data frame. Default is "`wt`".
#' @param ind_categorical Logical vector. (Optional) Indicator of which
#'.       `covariates` are categorical.
#'
#' @return A data frame of filtered matches
#' @export
#'
#' @examples
#' get_filtered_matches(
#'   example_trial_data,
#'   treatment = "Z",
#'   covariates = c("age", "stage4", "site", "ecog", "prevTrt", "diseaseFree")
#' )
get_filtered_matches <- function(data, treatment, covariates,
                                 id = NULL, weight = "wt", ind_categorical = NULL) {

  if (!(treatment %in% names(data))) {
    cli::cli_abort(c(
      "x" = "You input {treatment} for {.var treatment}, but it is not found in",
      " {.var data}."
    ))
  }

  missing_covs <- setdiff(covariates, names(data))
  if (length(missing_covs) > 0) {
    cli::cli_abort(c("x" = "Missing covariates in data: {missing_covs}"))
  }

  if (!is.null(ind_categorical)) {
    if (length(ind_categorical) != length(covariates)) {
      cli::cli_abort(c(
        "x" = "The length of {.var covariates} and {.var ind_categorical} must",
        " be equal."
      ))
    }
  }
  if (is.null(id)) {
    cli::cli_alert_info(c(
      "i" = "You did not provide a value for {.var id}, so we created a unique",
      " id called `id`"
    ))
    id <- "id"
    data$id <- seq_len(nrow(data))
  }

  dta_sorted <- data[order(data[[treatment]], decreasing=TRUE), ]

  dta_sorted$tID <- paste(dta_sorted[[treatment]], dta_sorted[[id]], sep="_")

  # Identify categorical variables
  if (is.null(ind_categorical)) {
    cov_data <- dta_sorted[, covariates, drop = FALSE]
    ind_categorical <- sapply(cov_data, is_categorical)
    cli::cli_alert_info(c(
      "i" = "You did not provide an indicator of categorical variables via",
      " {.var ind_categorical}, so we guessed that {covariates[ind_categorical]} are",
      " categorical."
    ))
  }

  cov_data <- dta_sorted[, covariates, drop = FALSE]
  if (any(ind_categorical == 1)) {
    cov_data[, which(ind_categorical == 1)] <- sapply(
      cov_data[, which(ind_categorical == 1), drop = FALSE],
      function(x) as.numeric(as.character(x)))
  }

  # Calculate Mahalanobis distance matrix
  dist <- get_mahalanobis_mat(X = cov_data)
  rownames(dist) <- dta_sorted$tID
  colnames(dist) <- dta_sorted$tID

  # Extract distance matrix between treated (rows) and control (columns)
  dist2 <- dist[substr(rownames(dist), 1, 1) == 1,
                substr(colnames(dist), 1, 1) == 0]

  # Transpose if more treated than controls
  if (dim(dist2)[1] > dim(dist2)[2]) {
    dist2 <- t(dist2)
  }

  # TODO find a better way to pass this option and don't supress all warnings
  options("optmatch_max_problem_size" = Inf)
  pm <- suppressWarnings(optmatch::pairmatch(dist2))

  match_df <- data.frame(id = names(pm), match_id = pm)
  match_df2 <- match_df[order(match_df$match_id, match_df$id), ]
  match_df <- match_df2[!is.na(match_df2$match_id), ]

  id_lookup <- stats::setNames(1:nrow(dta_sorted), dta_sorted$tID)

  match_groups <- split(match_df$id, match_df$match_id)

  result_list <- list()

  for (i in seq_along(match_groups)) {
    match_id <- names(match_groups)[i]
    pair <- match_groups[[i]]

    if (length(pair) != 2) {
      cli::cli_alert_warning(c(
        "Dropping matched pair {pair} for lack of match."
      ))
    } else {

      control_id <- pair[substr(pair, 1, 1) == "0"]
      treated_id <- pair[substr(pair, 1, 1) == "1"]

      if (length(control_id) == 1 || length(treated_id) == 1) {

        control_idx <- id_lookup[control_id]
        treated_idx <- id_lookup[treated_id]

        control_row <- dta_sorted[control_idx, , drop = FALSE]
        treated_row <- dta_sorted[treated_idx, , drop = FALSE]

        control_row$match_id <- treated_row$match_id <- match_id

        match_weight <- dist2[treated_id, control_id]
        control_row[[weight]] <- treated_row[[weight]] <- match_weight

        # Check if categorical variables match
        keep_match <- TRUE
        if (any(ind_categorical == 1)) {
          cat_vars <- covariates[ind_categorical == 1]
          for (var in cat_vars) {
            if (control_row[[var]] != treated_row[[var]]) {
              keep_match <- FALSE
              cli::cli_alert_warning(c(
                "Dropping pair {pair} because {var} did not exactly match ",
                "between groups."
              ))
            }
          }
        }
        if (keep_match) {
          result_list[[length(result_list) + 1]] <- rbind(control_row, treated_row)
        }
      }
    }
  }

  if (length(result_list) == 0) {
    cli::cli_alert_warning("No matches detected")
    return(data.frame())
  }

  matched_data <- do.call(rbind, result_list)

  matched_data$tID <- NULL

  return(matched_data)
}
