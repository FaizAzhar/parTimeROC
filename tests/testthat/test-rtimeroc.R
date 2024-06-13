test_that("Generating random data", {
  mod1 <- timeroc_obj(dist = "lognormal-weibull-PH")
  mod2 <- timeroc_obj(dist = "lognormal-weibull-copula", copula = "gaussian")
  NN <- 100
  x <- c(meanlog=1,sdlog=0.8)
  t <- c(shape=1.6,scale=1.2)
  ph <- 1
  cc <- -0.5
  expect_error(rtimeroc(n=NN, params.x=x, params.t=t, params.ph=ph), "Please supply a TimeROC object")
  expect_error(rtimeroc(mod1, params.x=x, params.t=t, params.ph=ph), "Please provide number of sample to be simulated")
  expect_error(rtimeroc(mod2, params.x=x, params.t=t, params.copula=cc), "Please provide number of sample to be simulated")
  expect_warning(rtimeroc(mod1, n=NN, params.t=t, params.ph=ph))
  expect_equal(all(is.na(rtimeroc(mod1, n=NN, params.x=x, params.ph=ph)$t)), TRUE)
  expect_error(rtimeroc(mod2, n = NN, params.x=x, params.t=t), "Please provide the parameter for copula")
  set.seed(12345)
  expect_equal(inherits(rtimeroc(mod1, n=NN, params.x=x, params.t=t, params.ph=ph), 'data.frame'), TRUE)
})
