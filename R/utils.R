is_categorical <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}

create_formula <- function(outcome, covariates) {
  covariates_string <- paste(covariates, collapse = " + ")
  formula_string <- paste(outcome, "~", covariates_string)
  formula_obj <- stats::as.formula(formula_string)
  return(formula_obj)
}
