library(testthat)
library(blblm)

#COEF
fit_parallelization <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, 4)

test_that("coefficient is a vector", {
  print(coef(fit_parallelization))
  expect_equal(length(coef(fit_parallelization)), 4)
}
)

#SIGMA
test_that("sigma is a vector",{
  print(sigma(fit_parallelization))
  expect_equal(length(sigma(fit_parallelization)), 1)
})

test_that("sigma with confidence is a vector",{
  print(sigma(fit_parallelization, confidence = TRUE))
  expect_equal(length(sigma(fit_parallelization, confidence = TRUE)), 3)
})

#CONFT
confint(fit_parallelization, c("wt", "hp"))


#PREDICT
test_that("predict with upper and lower bounds is a vector of length 6 (2 by 3)",{
  expect_equal(length(predict(fit_parallelization, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE)), 6)
})

test_that("predict is a vector of lengh 2",{
  expect_equal(length(predict(fit_parallelization, data.frame(wt = c(2.5, 3), hp = c(150, 170)))), 2)
})


