test_that("Fitting timeROC object", {
  mod1 <- timeroc_obj(dist = "lognormal-weibull-PH")
  mod2 <- timeroc_obj(dist = "lognormal-weibull-copula", copula = "gaussian")
  NN <- 100
  x <- c(meanlog=1,sdlog=0.8)
  t <- c(shape=1.6,scale=1.2)
  ph <- 1
  cc <- -0.5
  set.seed(12345)
  rr1 <- rtimeroc(mod1,n = NN, params.x = x, params.t = t, params.ph = ph)
  rr2 <- rtimeroc(mod2,n = NN, params.x = x, params.t = t, params.copula = cc)
  expect_equal(inherits(timeroc_fit(mod1, x = rr1$x, t = rr1$t, event = rr1$event), 'fitTROC'), TRUE)
  expect_equal(inherits(timeroc_fit(mod2, x = rr2$x, t = rr2$t, event = rr2$event), 'fitTROC'), TRUE)
  expect_equal(inherits(timeroc_fit(mod1, x = rr1$x, t = rr1$t, event = rr1$event, method = 'bayes'), 'fitTROC'), TRUE)
  expect_equal(inherits(timeroc_fit(mod2, x = rr2$x, t = rr2$t, event = rr2$event, method = 'bayes'), 'fitTROC'), TRUE)
  expect_error(timeroc_fit(x = rr1$x, t = rr1$t, event = rr1$event), "Please supply a TimeROC object")
  expect_error(timeroc_fit(mod1, t = rr1$t, event = rr1$event), "Please provide data for X, T and Event")
  expect_error(timeroc_fit(mod1, x = rr1$x, event = rr1$event), "Please provide data for X, T and Event")
  expect_error(timeroc_fit(mod1, x = rr1$x, t = rr1$t), "Please provide data for X, T and Event")
})
