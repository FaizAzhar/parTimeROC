test_that("Calculating rate change of ROC", {
  test <- timeroc_obj("normal-weibull-PH",
  params.x = c(mean=5, sd=0.8),
  params.t = c(shape=1.6, scale=5),
  params.ph = 1)
  expect_equal(inherits(rate_change(test, t = c(.1,.2)), 'list'), TRUE)
  test <- timeroc_obj("normal-weibull-copula", copula = "gaussian")
  rr <- rtimeroc(test, n=500,
                 params.x = c(mean=5, sd=0.8),
                 params.t = c(shape=1.6, scale=5),
                params.copula = -0.3)
  cc <- timeroc_fit(test, x = rr$x, t = rr$t, event = rr$event)
  expect_equal(inherits(rate_change(cc, t = c(1,2)), 'list'), TRUE)
  expect_error(rate_change(), "argument \"obj\" is missing, with no default")
  expect_error(rate_change(123), "Please supply a fitTROC or TimeROC object")
})
