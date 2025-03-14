#' Function to obtain decision tree
#'
#' @param data Data frame with averaged matches (obtained from [`get_averaged_matches()`])
#' @param covariates Character vector. Vector of covariate column names
#'        from `data` to use in match
#' @param survivial Logical. Indicator whether this is a survival analysis
#' @param l Numeric. Left
#' @param r Numeric. Right
#' @param weight Character. The name that should be given to the variable
#'        containing the weights in the output data frame. Default is "`wt`".
#' @param delta Character. The name that should be given to the delta variable
#'        in the output data frame. Default is "`delta`".
#' @param bins Numeric. The number of bins to use to categorize the weights.
#'        Default: `4`.
#'
#' @return A decision tree object
#' @export
#'
#' @examples
#' covariates <- c("age", "stage4", "site", "ecog", "prevTrt", "diseaseFree")
#' matches <- get_filtered_matches(
#'   example_trial_data,
#'   treatment = "Z",
#'   covariates = covariates
#' )
#' m <- get_averaged_matches(
#'   matches,
#'   treatment = "Z",
#'   outcome = "Y",
#'   covariates = covariates
#'   )
#' get_mtree(m, covariates)
get_mtree <- function(data, covariates, survivial = FALSE, l = -Inf, r = Inf,
                      weight = NULL, delta = "delta", bins = 4) {
  data$weight = ifelse(is.null(weight),
                       NA,
                       (bins + 1) -
                         as.numeric(Hmisc::cut2(data[[weight]], g = bins))
                       )
  if (survivial == FALSE) {
    if (dim(table(data[[delta]])) <= 3) {
      data[[delta]] <- as.factor(data[[delta]])
      }
    if (is.null(weight)) {
      fit <- partykit::ctree(create_formula(delta, covariates), data = data)
    } else {
      data[[weight]] = (bins + 1) -
        as.numeric(Hmisc::cut2(data[[weight]], g = bins))
      fit <- partykit::ctree(create_formula(delta, covariates), data = data,
                  weights = data[[weight]])
    }
  } else if (survivial) {
    data$time1[data$time1 == -Inf] <- l
    data$time2[data$time2 == Inf] <- r
    fit <- LTRCtrees::LTRCIT(
      create_formula("survival::Surv(time1, time2, type='interval2')",
                     covariates),
      data = data)
  }
  fit
}
