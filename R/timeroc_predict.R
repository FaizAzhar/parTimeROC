#' timeroc_predict
#'
#' @description Predict time-dependent ROC from fitted model.
#'
#' @param obj A 'fitTROC' or 'TimeROC' object.\cr
#' @param B An integer specifying bootstrap iteration. If B > 1, will also return confidence interval.\cr
#' @param cutoff A numeric specifying total cutoff point on ROC curve.\cr
#' @param t A numeric/vector specifying time point of interest. (Default: Time-to-event at 50th quantile points)\cr
#' @param newx A numeric/vector specifying biomarker of interest.\cr
#' @param type A string indicate type of analysis to run. (Default = 'standard')\cr
#' @param method A string specifying method of estimation. (Default = 'mle') \cr
#' @param definition A string indicating ROC definition to use. (Default = 'c/d') \cr
#' @param ci An integer 0 to 1 for confidence level.\cr
#' @param h An integer specifying small change of time (To compute density from S(t|x)) \cr
#' @param params.x A named vector for biomarker's parameter.\cr
#' @param params.t A named vector for time-to-event's parameter.\cr
#' @param copula A string indicating the type of copula to be used.\cr
#' @param params.copula An integer specifying the copula's parameter.\cr
#' @param params.ph An integer specifying the PH parameter.\cr
#' @param seed A numeric to pass in set.seed.\cr
#'
#' @examples
#' # PH model
#' test <- timeroc_obj('normal-weibull-PH')
#' set.seed(23456)
#' rr <- rtimeroc(obj = test, censor.rate = 0.1, n=500,
#'                params.t = c(shape=1, scale=5),
#'                params.x = c(mean=5, sd=1),
#'                params.ph=0.5)
#' cc <- timeroc_fit(x=rr$x, t=rr$t, event=rr$event, obj = test)
#' start.t <- Sys.time()
#' jj <- timeroc_predict(cc)
#' print(Sys.time()-start.t)
#' start.t <- Sys.time()
#' jj <- timeroc_predict(cc, t = quantile(rr$t,probs = seq(0,1,by=0.1)))
#' print(Sys.time()-start.t)
#'
#'
#' # Copula model
#' test <- timeroc_obj(dist = 'gompertz-gompertz-copula', copula='clayton90',
#' params.t = c(shape=3,rate=1),
#' params.x = c(shape=1,rate=2),
#' params.copula=-5)
#' set.seed(23456)
#' rr <- rtimeroc(obj = test, censor.rate = 0.2, n=500)
#' cc <- timeroc_fit(x=rr$x, t=rr$t, event=rr$event, obj = test)
#' start.t <- Sys.time()
#' jj <- timeroc_predict(cc)
#' print(Sys.time()-start.t)
#' start.t <- Sys.time()
#' jj <- timeroc_predict(cc, t = quantile(rr$t,probs = seq(0,1,by=0.1)))
#' print(Sys.time()-start.t)
#'
#' @returns A list of ROC dataframe for each time-to-event.
#' @export
timeroc_predict <- function(obj, t, newx, cutoff = 100, B = 1, type = 'standard',
                         params.x, params.t, copula,
                         method = 'mle', definition = 'c/d', seed,
                         params.copula, params.ph, ci=0.95, h=-0.0001){

  ## preprocessing of arguments
  se.x <- se.t <- se.c <- x.val <- name.dist.x <- name.dist.t <- dat <- NULL
  se.ph <- xval <- x.dist <- t.dist <- iscopula <- NULL

  if(inherits(obj, 'fitTROC')){
    args <- preproc(c(as.list(environment()), call = match.call()),
                  extract_from_fitTROC)
  }else if(inherits(obj, 'TimeROC')){
    args <- preproc(c(as.list(environment()), call = match.call()),
                    TimeROC_predict_process)
  }else{
    stop("Please supply a fitTROC or timeROC object")
  }

  list2env(args, environment())

  if(type == 'landmark'){
    res <- analysis.landmark(dat = dat, method = method, x.dist = x.dist, params.x = params.x,
                             ci = ci, t.dist = t.dist, params.t = params.t,
                             params.copula = params.copula, iscopula = iscopula, t = t,
                             copula = copula, se.x = se.x, se.t = se.t, se.c = se.c,
                             B = B, params.ph = params.ph, xval = xval, se.ph = se.ph,
                             definition = definition, h=h)
  } else if(type == 'standard'){
    res <- analysis.standard(iscopula = iscopula, params.x = params.x,
                             params.t = params.t, params.copula = params.copula,
                             params.ph = params.ph, se.x = se.x, se.t = se.t,
                             se.c = se.c, B = B, xval = xval,h=h,
                             t = t, x.dist = x.dist, t.dist = t.dist,
                             copula = copula, ci = ci, se.ph = se.ph, definition = definition)
  }

  class(res) <- "predictTROC"
  return(res)
}

# --------------------------------------
# function to run standard analysis.
analysis.standard <- function(iscopula, params.x, params.t, params.copula = NULL,h,
                              params.ph = NULL, se.x, se.t, se.c = NULL, se.ph = NULL,
                              B, xval, t, x.dist, t.dist, copula = NULL, ci, definition){

  if(iscopula == 1){
    pargs <- list(par.x = params.x, par.t = params.t,
                  par.cop = params.copula,
                  se.x = se.x, se.t = se.t, se.c = se.c)
    sim.par <- do.call(pars.boot,list(boot.val=B, pargs=pargs, iscopula=iscopula))
    ret <- do.call(timeroc,list(pars=sim.par, cutoff=xval, t=t, model=iscopula,
                                x.dist = x.dist, t.dist = t.dist, copula = copula,
                                ci = ci, definition = definition))
    for(i in 1:length(ret)){ret[[i]] <- cbind(ret[[i]], assoc = params.copula)}
  } else if (iscopula == 0){
    pargs <- list(par.x = params.x, par.t = params.t, par.ph = params.ph,
                  se.x = se.x, se.t = se.t, se.ph = se.ph)
    sim.par <- do.call(pars.boot,list(boot.val=B, pargs=pargs, iscopula=iscopula))
    ret <- do.call(timeroc,list(pars=sim.par, cutoff=xval, t=t, model=iscopula,
                                x.dist = x.dist, t.dist = t.dist,
                                ci = ci, definition = definition,h=h))
    for(i in 1:length(ret)){ret[[i]] <- cbind(ret[[i]], assoc = params.ph)}
  }
  return(ret)
}

# --------------------------------------
# function to run landmark analysis.
analysis.landmark <- function(dat, method, x.dist, params.x, ci, t.dist, params.t,
                              params.copula = NULL, iscopula, t, copula = NULL,
                              se.x, se.t, se.c = NULL, B, params.ph = NULL,
                              se.ph = NULL, xval, definition, h){
  ## ToDo: Future work for landmark analysis
  res <- c()
  fitted.x <- fit.x(x = dat$x, method = method, x.dist = x.dist,
                    init.param.x = params.x, ci = ci)

  if(iscopula == 1){
    fitted.t <- fit.t(x = dat$x, t = dat$t, event = dat$event,
                      method = method, t.dist = t.dist,
                      init.param.t = params.t,
                      iscopula = iscopula, ci = ci)

    for(ttime in t){
      subdat <- as.data.frame(dat)
      subdat[which(dat$t > ttime),3] <- 0

      fitted.copula <- fit.copula(x = subdat$x, t = subdat$t, event = subdat$event,
                                  method = method, res.x = fitted.x, res.t = fitted.t,
                                  x.dist = x.dist, t.dist = t.dist,
                                  init.param.copula = params.copula, copula = copula, ci = ci)

      pargs <- list(par.x = fitted.x$par, par.t = fitted.t$par,
                    par.cop = fitted.copula$par,
                    se.x = fitted.x$se, se.t = fitted.t$se, se.c = fitted.copula$se)

      sim.par <- do.call(pars.boot,list(boot.val=B, pargs=pargs, iscopula=iscopula))
      ret <- do.call(timeroc,list(pars=sim.par, cutoff=xval, t=ttime, model=iscopula,
                                  x.dist = x.dist, t.dist = t.dist, copula = copula,
                                  ci = ci, definition = definition))
      ret[[1]] <- cbind(ret[[1]], assoc = fitted.copula$par)
      res <- c(res,ret)
    }
  } else if (iscopula == 0){
    for(ttime in t){
      subdat <- dat
      subdat[which(dat$t > ttime),3] <- 0
      fitted.t <- fit.t(x = subdat$x, t = subdat$t, res.x = fitted.x, event = subdat$event,
                        method = method, x.dist = x.dist, t.dist = t.dist, init.param.t = params.t,
                        init.param.ph = params.ph, iscopula = iscopula, ci = ci)
      tname <- !(names(fitted.t$par) %in% c("beta"))
      pargs <- list(par.x = fitted.x$par, par.t = fitted.t$par[tname],
                    par.ph = fitted.t$par["beta"],
                    se.x = fitted.x$se, se.t = fitted.t$se[tname], se.ph = fitted.t$se['beta'])
      sim.par <- do.call(pars.boot,list(boot.val=B, pargs=pargs, iscopula=iscopula))
      ret <- do.call(timeroc,list(pars=sim.par, cutoff=xval, t=ttime, model=iscopula,
                                  x.dist = x.dist, t.dist = t.dist,
                                  ci = ci, definition = definition, h=h))
      ret[[1]] <- cbind(ret[[1]], assoc = fitted.t$par["beta"])
      res <- c(res,ret)
    }
  }
  return(res)
}

# --------------------------------------
# helper function to generate random parameters used in bootstrap process.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
pars.boot <- function(boot.val, pargs, iscopula){
  sim <- list()
  if(boot.val <= 1){
    sim$par.x <- matrix(rep(pargs$par.x,2),ncol=length(pargs$par.x), byrow=T)
    colnames(sim$par.x) <- names(pargs$par.x)
    sim$par.t <- matrix(rep(pargs$par.t,2),ncol=length(pargs$par.t), byrow=T)
    colnames(sim$par.t) <- names(pargs$par.t)
    if(iscopula==0){
      sim$par.b <- matrix(rep(pargs$par.ph,2),ncol=length(pargs$par.ph), byrow=T)
      colnames(sim$par.b) <- names(pargs$par.ph)
    } else {
      sim$par.c <- matrix(rep(pargs$par.cop,2),ncol=length(pargs$par.cop), byrow=T)
      colnames(sim$par.c) <- names(pargs$par.cop)
    }
    return(sim)
  } else{
      # sim$par.x <- rmvnorm(boot.val, pargs$par.x, diag(pargs$se.x))
      par.x <- matrix(NA, nrow = boot.val, ncol = length(pargs$se.x),
                      dimnames = list(NULL,names(pargs$se.x)))
      for(i in 1:ncol(par.x)){
        par.x[,i] <- rnorm(boot.val, mean = pargs$par.x[[i]], sd = pargs$se.x[[i]])
      }
      sim$par.x <- par.x

      par.t <- matrix(NA, nrow = boot.val, ncol = length(pargs$se.t),
                        dimnames = list(NULL,names(pargs$se.t)))
      for (i in 1:ncol(par.t)){
          par.t[,i] <- rnorm(boot.val, mean = pargs$par.t[[i]], sd = pargs$se.t[[i]])
      }
      sim$par.t <- par.t
      if(iscopula==0){
        # sim$par.t <- rmvnorm(boot.val, c(pargs$par.t,pargs$par.ph), diag(pargs$se.t))
        # sim$par.b <- sim$par.t[,'beta']
        # sim$par.t <- sim$par.t[,-'beta']
        sim$par.b <- rnorm(boot.val, mean = pargs$par.ph, sd = pargs$se.ph)
      } else {
        # sim$par.t <- rmvnorm(boot.val, c(pargs$par.t), diag(pargs$se.t))
        sim$par.c <- rnorm(boot.val, pargs$par.cop, pargs$se.c)
      }
      return(sim)
  }
}

# --------------------------------------
# helper function to compute sensitivity and specificity
timeroc <- function(model, pars, t, cutoff, x.dist, t.dist, copula=NULL, ci, definition, h){
  xfn <- x.dist$density
  tfn <- t.dist$density
  if(x.dist$name.abbr == "skewnormal") class(xfn) <- "snorm" # there is no sn::psn(..,lower.tail=F)
  if(t.dist$name.abbr == "skewnormal") class(tfn) <- "snorm"
  if(model == 1) { #iscopula = 1
    est <- do.call(roc.cop,list(x = cutoff, t = t, xfn = xfn, tfn = tfn,
                                cfn = copula, pars = pars, ci = ci))
  } else {
    est <- do.call(roc.ph,list(x = cutoff, t = t, xfn = xfn, tfn = tfn,
                               pars = pars, ci = ci, definition = definition, h=h))
  }

  return(est)
}

# --------------------------------------
# Function that will be automatically called within timeroc to compute sensitivity and specificity of a PH mdoel.
#' @importFrom cubature hcubature
#' @importFrom stats quantile na.omit
roc.ph <- function(x,t, xfn, tfn, pars, ci, definition, h){
  pobs <- (1-ci)/2
  res <- list()
  res.x <- matrix(NA, nrow = length(x), ncol = 7,
                  dimnames = list(NULL, c('est.sens', 'est.spec', 'low.sens',
                                          'low.spec', 'upp.sens', 'upp.spec','X')))
  res.boot <- matrix(NA, nrow = length(pars$par.b), ncol = 2)
  colnames(res.boot) <- c('sens','spec')

  Survt <- function(pars){
    st <- do.call(tfn[[3]], c(list(t_val,lower.tail = F),pars))
    if(inherits(tfn,"snorm")) st <- 1-st
    st
  }

  Survx <- function(pars){
    sx <- do.call(xfn[[3]], c(list(x_val,lower.tail = F),pars))
    if(inherits(tfn,"snorm")) sx <- 1-sx
    sx
  }

  Survt_h <- function(pars){
    st <- do.call(tfn[[3]], c(list(t_val+h,lower.tail = F),pars))
    if(inherits(tfn,"snorm")) st <- 1-st
    st
  }

  # For C/D time-dependent ROC
  # args[1] = x, args[2] = S_0(t), args[3] = beta
  integrate.f <- function(args){
    jj <- function(cc) {
      dx <- do.call(xfn[[2]], c(list(x = cc), args[5:length(args)]))
      return(args[2]^(exp(args[3] * cc)) * dx)
    }
    hcubature(jj, lowerLimit = args[1], upperLimit = Inf)$integral
  }

  # For I/D time-dependent ROC
  # args[4] = S_0(t+h)
  integrate.f.id <- function(args){
    jj <- function(cc) {
      dx <- do.call(xfn[[2]], c(list(x = cc), args[5:length(args)]))
      return((args[2]^(exp(args[3] * cc)) - (args[4])^(exp(args[3] * cc)))/(h) * dx)
    }
    hcubature(jj, lowerLimit = args[1], upperLimit = Inf)$integral
  }

  # # For C/D time-dependent ROC (Pepe definition)
  # # args[4] = S_0(t+h)
  # integrate.f.cd <- function(args){
  #   jj <- function(cc) {
  #     dx <- do.call(xfn[[2]], c(list(x = cc), args[5:length(args)]))
  #     return((args[2]^(exp(args[3] * cc)) - (args[4])^(exp(args[3] * cc))) * dx)
  #   }
  #   hcubature(jj, lowerLimit = args[1], upperLimit = Inf)$integral
  # }
  for(t_val in t){

    pt <- apply(pars$par.t,1,Survt)
    pt_h <- apply(pars$par.t,1,Survt_h)

    for(ii in seq_along(x)){
      x_val <- x[ii]
      px <- apply(pars$par.x,1,Survx)
      S.t <- apply(cbind(-Inf,pt,pars$par.b,pt_h,pars$par.x),1,integrate.f)
      s.ct <- apply(cbind(x_val,pt,pars$par.b,pt_h,pars$par.x),1,integrate.f)

      # specificity
      spec <- 1 - (s.ct/S.t)

      # sensitivity
      if(definition == "c/d"){
        # S.t.cd <- apply(cbind(x_val,pt,pars$par.b,pt_h,pars$par.x),1,integrate.f.cd)
        # S.t <- apply(cbind(-Inf,pt,pars$par.b,pt_h,pars$par.x),1,integrate.f.cd)
        sens <- (px-s.ct) / (1-S.t)
        # sens <- S.t.cd / S.t
      } else if(definition == "i/d"){
        f.t <- apply(cbind(-Inf,pt,pars$par.b,pt_h,pars$par.x),1,integrate.f.id)
        S.t.id <- apply(cbind(x_val,pt,pars$par.b,pt_h,pars$par.x),1,integrate.f.id)
        sens <- S.t.id/f.t
      }

      res.boot <- na.omit(cbind(sens, spec))
      est <- t(colMeans(res.boot))
      low <- t(apply(res.boot, 2, quantile, probs = pobs))
      upp <- t(apply(res.boot, 2, quantile, probs = 1-pobs))
      res.x[ii,] <- cbind(est, low, upp, x_val)
    }

    res[[as.character(t_val)]] <- res.x

  }

  return(res)
}

# --------------------------------------
# roc.cop
# Function that will be automatically called within timeroc to compute sensitivity and specificity of a copula mdoel.
#' @importFrom stats quantile na.omit
roc.cop <- function(x,t, xfn, tfn, cfn, pars, ci){
  pobs <- (1-ci)/2
  res <- list()
  res.x <- matrix(NA, nrow = length(x), ncol = 7,
                  dimnames = list(NULL, c('est.sens', 'est.spec', 'low.sens',
                                          'low.spec', 'upp.sens', 'upp.spec','X')))
  res.boot <- matrix(NA, nrow = length(pars$par.b), ncol = 2)
  colnames(res.boot) <- c('sens','spec')

  Ft <- function(pars){
    do.call(tfn[[3]], c(list(t_val),pars))
  }

  Fx <- function(pars){
    do.call(xfn[[3]], c(list(x_val),pars))
  }

  roc.formula <- function(args){
    F.uv <- cfn$density[[3]](u1 = args[1], u2 = args[2],
                             family = cfn$family, par = args[3])
    sens <- (args[2] - F.uv)/args[2]
    spec <- (args[1] - F.uv)/(1-args[2])
    data.frame(sens, spec)
  }

  for(t_val in t){
    v <- apply(pars$par.t,1,Ft)

    for(ii in seq_along(x)){
      x_val <- x[ii]
      u <- apply(pars$par.x,1,Fx)

      res.boot <- apply(cbind(u,v,pars$par.c),1,roc.formula)
      res.boot <- matrix(unlist(res.boot),ncol=2, nrow=length(pars$par.c),
                         byrow = TRUE)
      est <- t(na.omit(colMeans(res.boot)))
      low <- t(apply(res.boot, 2, quantile, probs = pobs, na.rm = T))
      upp <- t(apply(res.boot, 2, quantile, probs = 1-pobs, na.rm = T))

      res.x[ii,] <- cbind(est, low, upp, x_val)
    }

    res[[as.character(t_val)]] <- res.x

  }

  return(res)
}
