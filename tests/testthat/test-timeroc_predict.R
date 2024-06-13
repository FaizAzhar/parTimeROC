test_that("Predicting time-ROC", {
  mod1 <- timeroc_obj(dist = "lognormal-weibull-PH")
  NN <- 100
  x <- c(meanlog=1,sdlog=0.8)
  t <- c(shape=1.6,scale=1.2)
  ph <- 1
  mod2 <- timeroc_obj(dist = "lognormal-weibull-PH", params.x =x, params.t = t, params.ph = ph)
  set.seed(12345)
  rr1 <- rtimeroc(mod1,n = NN, params.x = x, params.t = t, params.ph = ph)
  ff <- timeroc_fit(mod1, x = rr1$x, t = rr1$t, event = rr1$event)
  expect_error(timeroc_predict(mod1))
  expect_equal(inherits(timeroc_predict(ff), 'predictTROC'), TRUE)
  expect_equal(inherits(timeroc_predict(mod2), 'predictTROC'), TRUE)
  expect_error(timeroc_predict(ff, t = 'abcd'), 'Non-numeric argument to mathematical function')
  expect_error(timeroc_predict(rr1), "Please supply a fitTROC or timeROC object")
  expect_equal(inherits(timeroc_predict(ff, newx = c(1,2,3)), 'predictTROC'), TRUE)
})
