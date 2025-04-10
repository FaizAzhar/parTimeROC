#' timeroc_fit
#'
#' @description Fit TimeROC using Maximum Likelihood Estimator.
#'
#' @param x A numeric vector of single biomarker or covariate.\cr
#' @param t A numeric vector of time-to-event.\cr
#' @param event A numeric vector of event status (0=dead, 1=alive).\cr
#' @param obj An initialized 'TimeROC' object.\cr
#' @param init.param.x Vector of starting value for biomarker parameter.\cr
#' @param init.param.t Vector of starting value for time-to-event parameter.\cr
#' @param init.param.copula An integer of starting value for copula parameter.\cr
#' @param init.param.ph An integer of starting value for association parameter.\cr
#' @param init_bayes Starting value when running the bayesian estimation.\cr
#' @param ci An integer 0 to 1 for confidence level.\cr
#' @param method A string specifying method of estimation. (Default = 'mle') \cr
#' @param weights Weights to handle Inverse Probability Censoring Weights. \cr
#' @param breakpoints Break points to specify intervals for Piecewise Hazard Model. \cr
#' @returns return a list of frequentist or bayesian estimator.
#' @examples
#' ## fitting copula model
#' test <- timeroc_obj(dist = 'gompertz-gompertz-copula', copula = "gumbel90")
#' set.seed(23456)
#' rr <- rtimeroc(obj = test, censor.rate = 0, n=500,
#'                params.t = c(shape=3,rate=1),
#'                params.x = c(shape=1,rate=2),
#'                params.copula=-5) # name of parameter must follow standard
#'
#' plot(t~x, rr)
#' start.t <- Sys.time()
#' cc <- timeroc_fit(rr$x, rr$t, rr$event, obj = test)
#' print(Sys.time()-start.t)
#'
#' ## fitting PH model
#' test <- timeroc_obj(dist = 'weibull-lognormal-PH')
#' set.seed(23456)
#' rr <- rtimeroc(obj = test, censor.rate = 0, n=100,
#'                params.t = c(meanlog=0, sdlog=1),
#'                params.x = c(shape=2, scale=1),
#'                params.ph=0.5) # name of parameter must follow standard
#' plot(t~x, rr)
#' start.t <- Sys.time()
#' cc <- timeroc_fit(rr$x, rr$t, rr$event, obj = test)
#' print(Sys.time()-start.t)
#'
#' @export
#' @importFrom stats qnorm glm confint.default
timeroc_fit <- function(obj, x, t, event, init.param.x = NULL, init.param.t= NULL,
                    init.param.copula = NULL, init.param.ph = NULL, ci = 0.95,
                    method = "mle", weights = NULL, breakpoints = NULL, init_bayes = NULL){

  if (missing(x) | missing(t) | missing(event)) stop("Please provide data for X, T and Event")

  ## preprocessing of arguments
  copula <- x.dist <- t.dist <- NULL
  args <- preproc(c(as.list(environment()), call = match.call()),
                  extract_from_TimeROC)
  list2env(args, environment())

  if (is.null(init.param.ph) & !is.list(copula)) init.param.ph <- survival::coxph(survival::Surv(t,event) ~ x - 1,
                                                            data = data.frame(x=x,t=t,event=event))$coef[[1]]

  if (is.null(init.param.x)) init.param.x <- x.dist$init(x) else as.list(init.param.x)
  if (is.null(init.param.t) & !is.list(copula)) {
      if(t.dist$name.abbr != 'pch'){
        init.param.t <- t.dist$init(t,x,event)
      } else if(t.dist$name.abbr == 'pch'){
        init.param.t <- t.dist$init(breakpoints)$rates
        t.dist$lower <- rep(0,length(init.param.t))
        t.dist$upper <- rep(Inf, length(init.param.t))
      }
    }
  else if (is.null(init.param.t) & is.list(copula)){init.param.t <- t.dist$init(t)}
  else as.list(init.param.t)

  ## Fit x
  res.x <- fit.x(x, method, x.dist, init.param.x, ci)

  ## Fit t
  if(is.null(weights)) weights <- rep(1,length(t))
  if(t.dist$name.abbr != 'pch'){
    res.t <- fit.t(x, t, res.x, event, method, x.dist, t.dist,
                   init.param.t, init.param.ph, is.list(copula), ci,
                   weights, breakpoints, init_bayes)
  }else if(t.dist$name.abbr == "pch"){
    if(is.null(breakpoints)) stop("Please specify time cutoff")
    res.t <- fit.tpch(data.frame(x,t,event), breakpoints, ci)
  }

  ##Fit Copula
  if (is.list(copula)){
    res.c <- fit.copula(x, t, res.x, event, method, res.t, x.dist, t.dist, init.param.copula, copula, ci)

    out <- list(name = obj$dist, name.copula = obj$copula$cop.abbr,
                dat = data.frame(x=x,t=t,event=event), x = res.x, t = res.t,
                copula = res.c, conf = ci)
    class(out) <- "fitTROC"
    return(out)
  }

  out <- list(name = obj$dist, dat = data.frame(x=x,t=t,event=event),
              x = res.x, t = res.t, conf = ci)
  class(out) <- "fitTROC"
  return(out)
}

# Routine to fit copula
fit.copula <- function(x, t, res.x, event, method, res.t, x.dist, t.dist, init.param.copula, copula, ci){
  Call3 <- match.call(expand.dots = TRUE)
  res.c <- list()
  xargs <- as.list(res.x$par)
  targs <- as.list(res.t$par)

  px <- as.call(c(list(x.dist$density[[3]], x), xargs))
  pt <- as.call(c(list(t.dist$density[[3]], t), targs))

  u <- eval(px)
  v <- eval(pt)

  if (is.null(init.param.copula)) init.param.copula <- copula$init(u,v)

  dx <- as.call(c(list(x.dist$density[[2]], x = x, log=TRUE), xargs))
  dt <- as.call(c(list(t.dist$density[[2]], x = t, log=TRUE), targs))
  log.dx <- eval(dx)
  log.dt <- eval(dt)

  ##Fit Copula
  if (method == "mle" | method == "bayes"){
    ll.c <- function(parms, ...){
      likeli <- event*(log(copula$density[[2]](u, v,
                                               family = copula$family,
                                               par = parms))+log.dx + log.dt) +
        (1-event)*(log(1 - copula$density[[4]](u,v,
                                               family = copula$family,
                                               par = parms)) + log.dx)
      -sum(likeli)
    }

    Call3[[1L]] <- quote(stats::optim)
    Call3$par <- init.param.copula
    Call3$fn <- ll.c
    Call3$hessian <- TRUE
    Call3$lower <- copula$lower
    Call3$upper <- copula$upper
    Call3$method <- "L-BFGS-B"
    res.c <- eval.parent(Call3)
    # if(res.c$convergence > 0L) stop("copula optimization failed. Consider method = `bayes`")
    res.c$cov <- .hess_to_cov(res.c$hessian)
    se <- sqrt(diag(res.c$cov))
    res.c$ubound <- res.c$par + qnorm(1 - (1 - ci)/2)*se
    res.c$lbound <- res.c$par - qnorm(1 - (1 - ci)/2)*se
    res.c$se <- se
    names(res.c$par) <- 'theta'
    res.c$aic <- 2 * length(res.c$par) + 2 * res.c$value

  } else{
    res.c <- NA
  }
  return(res.c)
}

# Routine to fit time-to-event
fit.t <- function(x, t, res.x, event, method, x.dist, t.dist,
                  init.param.t, init.param.ph = NULL, iscopula, ci,
                  weights = NULL, breakpoints = NULL, init_bayes = NULL){

    Call2 <- match.call(expand.dots = TRUE)
    res.t <- list()
    xargs <- as.list(res.x$par)

    dx <- as.call(c(list(x.dist$density[[2]], x), xargs))
    T.end <- max(t)
    u <- eval(dx)
    u <- log(u)

    # nx <- length(res.x$par)

    ## Fit t
    if (method == "mle"){
      Call2[[1L]] <- quote(stats::optim)

      if (!iscopula){
        ll.t <- function(parms, ...){
          # dx <- as.call(c(list(x.dist$density[[2]], x), parms[1:nx]))
          # u <- eval(dx)

          likeli <- weights*(event * (log(t.dist$hazard(t,parms)) + x * parms[length(parms)]) -
            exp(x * parms[length(parms)]) * t.dist$cum.hazard(t, parms))
          -sum(likeli)
        }
      } else {
        ll.t <- function(parms, ...){
          names(parms) <- names(t.dist$pars)
          d <- as.call(c(list(t.dist$density[[2]],x=t,log=TRUE), as.list(parms)))
          -sum(eval(d))
        }
      }
      Call2$fn <- ll.t
      if (!iscopula) {
        init.param.t <- c(init.param.t, beta = init.param.ph)
        Call2$par <- unlist(init.param.t)
      } else {Call2$par <- unlist(init.param.t)}

      Call2$hessian <- TRUE
      Call2$lower <- c(t.dist$lower, -Inf)
      Call2$upper <- c(t.dist$upper, Inf)
      Call2$method <- "L-BFGS-B"

      res.t <- eval.parent(Call2)
      # if(res.t$convergence > 0L) stop("Time-to-event optimization failed. Consider method = `bayes`")
      if (!iscopula) {
          names(res.t$par) <- c(t.dist$pars,'beta')
      } else {names(res.t$par) <- t.dist$pars}
      res.t$cov <- .hess_to_cov(res.t$hessian)
      se <- sqrt(diag(res.t$cov))
      res.t$ubound <- res.t$par + qnorm(1 - (1 - ci)/2)*se
      res.t$lbound <- res.t$par - qnorm(1 - (1 - ci)/2)*se
      res.t$se <- se
      res.t$aic <- 2 * length(res.t$par) + 2 * res.t$value
    } else if (method == "bayes" & !iscopula){
        p <- 1
        baseline_t <- set_baseline(t.dist$name.abbr)
        stan_data <- list(time = t, event = event,
                          X = x, n = length(x),
                          p = p, baseline = baseline_t)
        if(is.null(init_bayes)) fitted.t <- rstan::sampling(stanmodels$ph, data = stan_data)
        else {
          init_f1 <- function(){init_bayes}
          fitted.t <- rstan::sampling(stanmodels$ph, data = stan_data, init = init_f1)
        }
        summary.t <- rstan::summary(fitted.t, probs = c((1-ci)/2, 1-(1-ci)/2))
        NN <- nrow(summary.t[["summary"]])
        res.t$par <- summary.t[["summary"]][-NN,1]
        res.t$value <- summary.t[["summary"]][NN,1]
        res.t$ubound <- summary.t[["summary"]][-NN,5]
        res.t$lbound <- summary.t[["summary"]][-NN,4]
        res.t$se <- summary.t[["summary"]][-NN,3]
        res.t$n_eff <- summary.t[["summary"]][-NN,6]
        res.t$Rhat <- summary.t[["summary"]][-NN,7]
        names(res.t$ubound) <- names(res.t$lbound) <- names(res.t$se) <- names(res.t$n_eff) <- names(res.t$Rhat) <- names(res.t$par) <- c("beta",t.dist$pars)
        res.t$aic <- NA
    } else if (method == "bayes" & iscopula){
        baseline_t <- set_baseline(t.dist$name.abbr)
        stan_data <- list(X = t, n = length(t),
                          baseline = baseline_t)
        fitted.t <- rstan::sampling(stanmodels$marginal, data = stan_data)
        summary.t <- rstan::summary(fitted.t, probs = c((1-ci)/2, 1-(1-ci)/2))
        NN <- nrow(summary.t[["summary"]])
        res.t$par <- summary.t[["summary"]][-NN,1]
        res.t$value <- summary.t[["summary"]][NN,1]
        res.t$ubound <- summary.t[["summary"]][-NN,5]
        res.t$lbound <- summary.t[["summary"]][-NN,4]
        res.t$se <- summary.t[["summary"]][-NN,3]
        res.t$n_eff <- summary.t[["summary"]][-NN,6]
        res.t$Rhat <- summary.t[["summary"]][-NN,7]
        names(res.t$ubound) <- names(res.t$lbound) <- names(res.t$se) <- names(res.t$n_eff) <- names(res.t$Rhat) <- names(res.t$par) <- t.dist$pars
        res.t$aic <- NA
    }

  return(res.t)
}

# Routine to fit biomarker
fit.x <- function(x, method, x.dist, init.param.x, ci){

  Call <- match.call(expand.dots = TRUE)
  res.x <- list()

  ## Fit x
  if (method == "mle"){
    ll.x <- function(parms, ...){
      names(parms) <- names(x.dist$pars)
      d <- as.call(c(list(x.dist$density[[2]],x=x,log=TRUE), as.list(parms)))
      -sum(eval(d))
    }

    Call[[1L]] <- quote(stats::optim)
    Call$par <- unlist(init.param.x)
    Call$fn <- ll.x
    Call$hessian <- TRUE
    Call$lower <- x.dist$lower
    Call$upper <- x.dist$upper
    if (length(unlist(init.param.x)) > 1L) Call$method <- "L-BFGS-B"
    else Call$method <- "Brent"

    res.x <- eval.parent(Call)
    if(res.x$convergence > 0L) stop("biomarker optimization failed. Consider method = `bayes`")
    res.x$cov <- .hess_to_cov(res.x$hessian)
    se <- sqrt(diag(res.x$cov))
    res.x$ubound <- res.x$par + qnorm(1 - (1 - ci)/2)*se
    res.x$lbound <- res.x$par - qnorm(1 - (1 - ci)/2)*se
    res.x$se <- se
    res.x$aic <- 2 * length(res.x$par) + 2 * res.x$value
  } else if (method == "bayes"){
    baseline_x <- set_baseline(x.dist$name.abbr)
    stan_data <- list(X = x, n = length(x),
                      baseline = baseline_x)
    fitted.x <- rstan::sampling(stanmodels$marginal, data = stan_data)
    summary.x <- rstan::summary(fitted.x, probs = c((1-ci)/2, 1-(1-ci)/2))
    NN <- nrow(summary.x[["summary"]])
    res.x$par <- summary.x[["summary"]][-NN,1]
    res.x$value <- summary.x[["summary"]][NN,1]
    res.x$ubound <- summary.x[["summary"]][-NN,5]
    res.x$lbound <- summary.x[["summary"]][-NN,4]
    res.x$se <- summary.x[["summary"]][-NN,3]
    res.x$n_eff <- summary.x[["summary"]][-NN,6]
    res.x$Rhat <- summary.x[["summary"]][-NN,7]
    names(res.x$ubound) <- names(res.x$lbound) <- names(res.x$se) <- names(res.x$n_eff) <- names(res.x$Rhat) <- names(res.x$par) <- x.dist$pars
    res.x$aic <- NA
  }
  return(res.x)
}

# Routine to fit Piecewise Hazard/Exponential Models
setup.pch <- function(df, breakpoints){
  breakpoints <- c(0,breakpoints,Inf)
  pos <- findInterval(df$t, vec = breakpoints)
  exposure_time <- numeric(0)
  # Create event matrix
  mat_event <- matrix(0, nrow = nrow(df), ncol = (length(breakpoints)-1))
  for(i in seq_along(df$t)){
    mat_event[i,pos[i]] <- 1
    exposure_time[i] <- df$t[i] - breakpoints[pos[i]]
  }

  # Create time exposure matrix
  mat_exposure <- matrix(1:(length(breakpoints)-1), nrow = nrow(df),
                         ncol = (length(breakpoints)-1), byrow = TRUE)
  full_interval <- diff(breakpoints)
  for(i in seq_along(df$t)){
    mat_exposure[i,] <- mat_exposure[i,] - pos[i]
  }

  mat_exposure[mat_exposure > 0] <- 0
  mat_exposure[mat_exposure < 0] <- full_interval[col(mat_exposure)][mat_exposure < 0]
  for(i in seq_along(df$t)){
    mat_exposure[i,pos[i]] <- exposure_time[i]
  }

  return(list(mat_event = mat_event, mat_exposure = mat_exposure))
}

fit.tpch <- function(df, breakpoints, ci){
  breakpoints <- breakpoints[breakpoints != 0 & !is.infinite(breakpoints)]
  breakpoints <- sort(breakpoints)
  labels_break <- paste("(",c(0,breakpoints),",",c(breakpoints,Inf),"]")
  params_name <- c()

  for (i in seq_along(labels_break)){params_name <- append(params_name, paste0('lambda',i))}

  aug_dat <- survival::survSplit(Surv(t,event)~x,data = df,
                                 cut = breakpoints, start = "start", episode = "interval")

  aug_dat$exposure <- aug_dat$t - aug_dat$start
  aug_dat$interval <- factor(aug_dat$interval, labels = labels_break)

  # Fitting GLM
  fit <- glm(event ~ interval + x + offset(log(exposure)), data = aug_dat, family = "poisson")
  b <- coef(fit)
  conf_b <- confint.default(fit, level = ci)
  res.t <- list()
  res.t$par <- c(exp(b[1] + c(0,b[2:length(labels_break)])), b[(length(labels_break)+1):length(coef(fit))])
  names(res.t$par) <- append(params_name, 'beta')
  res.t$value <- sn::logLik(fit)[1]
  res.t$counts <- fit$iter
  res.t$convergence <- fit$converged
  res.t$hessian <- solve(sn::vcov(fit))
  res.t$cov <- .hess_to_cov(res.t$hessian)
  se <- sqrt(diag(res.t$cov))
  res.t$ubound <- c(exp(conf_b[1,2] + c(0,conf_b[2:length(labels_break),2])), conf_b[(length(labels_break)+1):length(coef(fit)),2])
  res.t$lbound <- c(exp(conf_b[1,1] + c(0,conf_b[2:length(labels_break),1])), conf_b[(length(labels_break)+1):length(coef(fit)),1])
  res.t$se <- se
  names(res.t$se) <- names(res.t$ubound) <- names(res.t$lbound) <- names(res.t$par) <- append(params_name, 'beta')
  res.t$aic <- fit$aic
  res.t$breakpoints <- c(0,breakpoints,Inf)
  return(res.t)
}

# helper function to safely convert a Hessian matrix to covariance matrix
#' @importFrom Matrix nearPD
.hess_to_cov <- function(hessian, tol.solve = 1e-16, tol.evalues = 1e-5, ...) {
  if(is.null(tol.solve)) tol.solve <- .Machine$double.eps
  if(is.null(tol.evalues)) tol.evalues <- 1e-5
  # use solve(.) over chol2inv(chol(.)) to get an inverse even if not PD
  # less efficient but more stable
  inv_hessian <- solve(hessian, tol = tol.solve)
  if (any(is.infinite(inv_hessian)))
    stop("Inverse Hessian has infinite values.  This might indicate that the model is too complex to be identifiable from the data")
  evalues <- eigen(inv_hessian, symmetric = TRUE, only.values = TRUE)$values
  if (min(evalues) < -tol.evalues)
    warning(sprintf(
      "Hessian not positive definite: smallest eigenvalue is %.1e (threshold: %.1e). This might indicate that the optimization did not converge to the maximum likelihood, so that the results are invalid. Continuing with the nearest positive definite approximation of the covariance matrix.",
      min(evalues), -tol.evalues
    ))
  # make sure we return a plain positive definite symmetric matrix
  as.matrix(Matrix::nearPD(inv_hessian, ensureSymmetry = TRUE, ...)$mat)
}
