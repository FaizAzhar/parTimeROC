test_that("Loading timeROC object", {
  expect_equal(inherits(timeroc_obj(dist = "lognormal-weibull-PH"), "TimeROC"), TRUE)
  expect_equal(inherits(timeroc_obj(dist = "lognormal-weibull-copula",
                                    copula = "gaussian"), "TimeROC"), TRUE)
  expect_error(timeroc_obj(), "Please provide the distributions and models to be used.
                          EG: normal-weibull-PH / normal-weibull-copula")
  expect_error(timeroc_obj(dist = "lognormal-weibull-copula"), "Please provide a listed copula")
})
