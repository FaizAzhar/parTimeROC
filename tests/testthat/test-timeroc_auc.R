test_that("Calculating time-AUC", {
  x <- c(meanlog=1,sdlog=0.8)
  t <- c(shape=1.6,scale=1.2)
  ph <- 1
  mod2 <- timeroc_obj(dist = "lognormal-weibull-PH", params.x =x, params.t = t, params.ph = ph)
  ff <- timeroc_predict(mod2)
  expect_error(timeroc_auc(mod2), "Please provide a predictTROC object")
})
