#' Calculate Effect Statistics
#'
#' @param data Data frame.
#' @param treatment Character. Name of the treatment column in `data`
#' @param outcome Character. Name of the outcome column in `data`
#' @param covariates Character vector. Vector of covariate column names
#'        from `data` to use in match
#' @param tree Decision tree obtained from [`get_mtree()`].
#' @param cens_var Character. Name of a censoring column in `data`. If none
#'        exists, leave `NULL`. Default: `NULL`
#'
#' @return Data frame with effect estimates
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
#' tree <- get_mtree(m, covariates)
#' get_mtree_stats(example_trial_data,
#'   treatment = "Z",
#'   outcome = "Y",
#'   covariates = covariates,
#'   tree = tree)
get_mtree_stats = function(data,
                           treatment,
                           outcome,
                           covariates,
                           tree,
                           cens_var = NULL) {
  ind.cat <- sapply(data[, covariates, drop = FALSE], is_categorical)

  # Make sure variable types are appropriate
  FAC <- covariates[which(ind.cat == 1)]
  if (length(FAC) > 0) {
    for (i in 1:length(FAC)) {
      data[, names(data) == FAC[i]] = as.factor(data[, names(data) == FAC[i]])
    }
  }
  NUM <- covariates[which(ind.cat == 0)]
  if (length(NUM) > 0) {
    for (i in 1:length(NUM)) {
      data[, names(data) == NUM[i]] = as.numeric(data[, names(data) == NUM[i]])
    }
  }

  # Get subgroup label for original observations
  data$predNode = stats::predict(tree, newdata = data, type = "node")

  #Look at balance within discovered subgroups
  data$Z <- data[, names(data) == treatment]
  #tabs <- halfmoon::tidy_smd(data, .vars = covariates, .group = treatment)
  data$Y = data[, names(data) == outcome]
  if (length(unique(data$Y)) <= 2 & is.null(cens_var)) {
    est = by(data, data$predNode, function(x)
      get_binary_est(x))
  } else if (length(unique(data$Y)) > 2 & is.null(cens_var)) {
    est = by(data, data$predNode, function(x)
      get_continuous_est(x))
  } else {
    data$event = data[[cens_var]]
    L = -max(data$Y)
    U = max(data$Y)
    est = by(data, data$predNode, function(x)
      get_survival_est(x, L, U))
  }
  output.all = list(
    #tabs = tabs,
    est = est)
  return(output.all)
}



get_binary_est = function(e1) {
  e1$Y <- as.numeric(as.character(e1$Y))
  p1 <- mean(e1$Y[e1$Z == 1])
  p0 <- mean(e1$Y[e1$Z == 0])
  ee1 <- matrix(table(-e1$Z,-e1$Y), 2, 2, byrow = FALSE)
  r2 = RI2by2::Perm.CI.RLH(ee1, 0.05, total_tests = 1000)
  data.frame(
    x1 = ee1[1, 1],
    n1 = sum(ee1[1,]),
    p1 = p1,
    x0 = ee1[2, 1],
    n0 = sum(ee1[2,]),
    p0 = p0,
    diff = r2$tau.hat,
    lower = r2$RLH[1],
    upper = r2$RLH[2]
  )
}


get_continuous_est = function(e1) {
  j = stats::wilcox.test(e1$Y[e1$Z == 1], e1$Y[e1$Z == 0], conf.int = TRUE)
  data.frame(
    diff = j$estimate,
    lower = j$conf.int[1],
    upper = j$conf.int[2]
  )
}



get_survival_est = function(e1, l, u, l2 = -10, u2 = 3) {
  est <- seq(l, u, (u - l) / 1000)
  pval <- double()
  for (i in 1:length(est)) {
    dta2 = e1
    dta2$Y[dta2$Z == 0] = dta2$Y[dta2$Z == 0] + est[i]
    j2 <-
      survival::survdiff(survival::Surv(Y, event) ~ Z,
                         data = dta2,
                         rho = 0) #logrank
    pval <- c(pval, j2$pvalue)
  }
  add <- data.frame(diff = est[which.max(pval)],
                    lower = min(est[pval > 0.05]),
                    upper = max(est[pval > 0.05]))
  est2 <- exp(seq(l2, u2, (u2 - l2) / 1000))
  pval2 <- double()
  for (i in 1:length(est2)) {
    dta2 <- e1
    dta2$Y[dta2$Z == 0] <- dta2$Y[dta2$Z == 0] * est2[i]
    j2 <-
      survival::survdiff(survival::Surv(Y, event) ~ Z,
                         data = dta2,
                         rho = 0) #logrank
    pval2 = c(pval2, j2$pvalue)
  }
  mult <- data.frame(diff = est2[which.max(pval2)],
                     lower = min(est2[pval2 > 0.05]),
                     upper = max(est2[pval2 > 0.05]))
  output.all = list(add = add, mult = mult)
  return(output.all)
}
