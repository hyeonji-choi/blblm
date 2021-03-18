#' @title blblogreg
#' @import purrr
#' @import stats
#' @import parallel
#' @import future
#' @import furrr
#' @import tidyverse
#' @importFrom magrittr %>%
#' @importFrom utils capture.output 
#' @aliases blblm-package
#' @aliases _PACKAGE
#' @details this is a package that calculate logistic regression with a little bag of bootstraps.
#' Logistic Regression with Little Bag of Bootstraps
"_PACKAGE"
#' NULL

utils::globalVariables(c("."))
## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R

#' @title blblogreg
#' @param formula  the formula of choice - logistic linear regression with weight in this case.
#' @param data data of choice
#' @param B B a posive integer. It represents a number of bootstrap samples.
#' @param ncpu Parallel is set as TRUE in default which means the algorithm is implementing multiple CPUs, 4. The user can select FALSE to implement only one CPU.
#' @param m a positive number, the number of items to choose from.
#' @export
blblogreg <- function(formula, data, m = 10, B = 5000, ncpu = 1) {
  if (!ncpu){
    suppressWarnings(plan(multiprocess, workers = ncpu))
    options(future.rng.onMisuse = "ignore")
    
    data_list <- split_data(data, m)
    estimates <- future_map(
      data_list,
      ~ logreg_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
    res <- future_map(
      data_list,
      ~ list(estimates = estimates, formula = formula)
    )
    class(res) <- "blblogreg"
    invisible(res)
  }
  else{
    data_list <- split_data(data, m)
    estimates <- map(
      data_list,
      ~ logreg_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblogreg"
    invisible(res)
  }
}

#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' @title logistic regression of each subsample
#' @param formula  the formula of choice - logistic regression with weight in this case.
#' @param data data of choice
#' @param B B a posive integer. It represents a number of bootstrap samples.
#' @param n a positive number, the number of samples.
#' @return linear regression of each subsample with bootstrap
#' compute the estimates
logreg_each_subsample <- function(formula, data, n, B) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(data))))
  replicate(B, logreg1(data, formula, freqs), simplify = FALSE)
}

#' @title logistic regression coefficient and sigma
#' @param data data
#' @param formula formula
#' @param freqs frequency
#' @return logistic regression
#' compute the regression estimates for a blb dataset
logreg1 <- function(data, formula, freqs) {
  fit <- glm(formula, family = "binomial", data, weights = freqs)
  list(coef = blbcoef(fit), odds = blbodd(fit))
}

#' @title coefficient of blb
#' @param fit fit
#' @return coefficient of the fit
#' compute the coefficients from fit
blbcoef <- function(fit) {
  coef(fit)
}

blbodd <- function(fit){
  exp(coef(fit))
}


#' @title printblblogreg
#' @method print blblogreg
#' @param x x values
#' @param ... unused variables
#' @export
print.blblogreg <- function(x, ...) {
  cat("blblogreg model:", capture.output(x$formula))
  cat("\n")
}

#' @export
#' @title sigmablblogreg
#' @method sigma blblogreg
#' @param object object
#' @param confidence whether or not the result shows upper and lower limits at 0.95 confidence level. Default is FALSE.
#' @param level this is a confidence level, default at 0.96=5.
#' @param ncpu positive interger, number of CPU. Default to be 1 which represents one CPU. 
#' @param ... unused variables
sigma.blblogreg <- function(object, confidence = FALSE, level = 0.95, ncpu = 1, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2)), ncpu) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @title coefblblogreg
#' @method coef blblogreg
#' @param object object
#' @param ncpu positive interger, number of CPU. Default to be 1 which represents one CPU. 
#' @param ... unused variables
coef.blblogreg <- function(object, ncpu = 1, ...) { #missing the parallel input
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans(), ncpu)
}


#' @export
#' @title confintblblogreg
#' @method confint blblogreg
#' @param object object
#' @param parm parameter. default is NULL
#' @param level this is a confidence level, default at 0.96=5.
#' @param ncpu positive interger, number of CPU. Default to be 1 which represents one CPU. 
#' @param ... unused variables
confint.blblogreg <- function(object, parm = NULL, level = 0.95, ncpu = 1,...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)), ncpu)
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @export
#' @title predictblblogreg
#' @method predict blblogreg
#' @param object object
#' @param confidence whether or not the result shows upper and lower limits at 0.95 confidence level. Default is FALSE.
#' @param level this is a confidence level, default at 0.96=5.
#' @param new_data new data that is being created.
#' @param ncpu positive interger, number of CPU. Default to be 1 which represents one CPU. 
#' @param ... unused variables
predict.blblogreg <- function(object, new_data, confidence = FALSE, level = 0.95, ncpu = 1, ...) {
  est <- object$estimatesk
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t(), ncpu)
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans(), ncpu)
  }
}
  
  
