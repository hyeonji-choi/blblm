#' @title blblm
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
#' @details this is a package the calculate linear regression with a little bag of bootstraps.
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"
#' NULL

utils::globalVariables(c("."))
## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R


#' @title blblm
#' @param formula  the formula of choice - linear regression with weight in this case.
#' @param data data of choice
#' @param B B a posive integer. It represents a number of bootstrap samples.
#' @param ncpu set as 1 in default which means the algorithm is implementing 1 CPU. The user can input any number of CPUs they would like to use.
#' @param m a positive number, the number of items to choose from.
#' @return linear regression model
#' @export
blblm <- function(formula, data, m = 10, B = 5000, ncpu = 1) {
  if (!ncpu){
    suppressWarnings(plan(multiprocess, workers = ncpu)) #make it an option
    options(future.rng.onMisuse = "ignore")
    data_list <- split_data(data, m)
    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
    res <- future_map(
      data_list,
      ~ list(estimates = estimates, formula = formula)
    )
    class(res) <- "blblm"
    invisible(res)
  }
  else{
    data_list <- split_data(data, m)
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblm"
    invisible(res)
  }
}

#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' @title lmeachsubsample
#' @param formula  the formula of choice - linear regression with weight in this case.
#' @param data data of choice
#' @param B B a posive integer. It represents a number of bootstrap samples.
#' @param n a positive number, the number of samples.
#' @return linear regression of each subsample
#' compute the estimates
lm_each_subsample <- function(formula, data, n, B) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  replicate(B, lm1(X, y, n), simplify = FALSE)
}

#' @title lm1
#' @param X x values
#' @param y y values
#' @param n number of samples
#' @return linear regression with weight, more specifically coefficiemt and sigma. 
#' compute the regression estimates for a blb dataset
lm1 <- function(X, y, n) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit <- lm.wfit(X, y, freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' @title coefficient of blb
#' @param fit fit
#' @return coefficient of the fit
#' compute the coefficients from fit
blbcoef <- function(fit) {
  coef(fit)
}

#' @title sigma of blb
#' @param fit fit
#' @return sigma of the fit
#' compute sigma from fit
blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @title printblblm
#' @method print blblm
#' @param x x values
#' @param ... unused variables
#' @export
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}



#' @title sigmablblm
#' @method sigma blblm
#' @param object object
#' @param confidence whether or not the result shows upper and lower limits at 0.95 confidence level. Default is FALSE.
#' @param level this is a confidence level, default at 0.96=5.
#' @param ncpu positive interger, number of CPU. Default to be 1 which represents one CPU. 
#' @param ... unused variables
#' @export
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @title coefblblm
#' @method coef blblm
#' @param object object
#' @param ncpu positive interger, number of CPU. Default to be 1 which represents one CPU. 
#' @param ... unused variables
#' @export
coef.blblm <- function(object, ...) { #missing the parallel input
    est <- object$estimates
    map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
  }

#' @title confintblblm
#' @method confint blblm
#' @param object object
#' @param parm parameter. default is NULL
#' @param level this is a confidence level, default at 0.96=5.
#' @param ncpu positive interger, number of CPU. Default to be 1 which represents one CPU. 
#' @param ... unused variables
#' @export
confint.blblm <- function(object, parm = NULL, level = 0.95,...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @title predictblblm
#' @method predict blblm
#' @param object object
#' @param confidence whether or not the result shows upper and lower limits at 0.95 confidence level. Default is FALSE.
#' @param level this is a confidence level, default at 0.96=5.
#' @param new_data new data that is being created.
#' @param ncpu positive interger, number of CPU. Default to be 1 which represents one CPU. 
#' @param ... unused variables
#' @export
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t(), ncpu)
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}

mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
    (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
    }

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}

