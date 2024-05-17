#' timeroc_predict
#'
#' @description Predict time-dependent ROC from fitted model.
#'
#' @param obj A 'fitTROC' object.\cr
#' @param B An integer specifying bootstrap iteration. If B > 1, will also return confidence interval.\cr
#' @param cutoff A numeric specifying total cutoff point on ROC curve.\cr
#' @param t A numeric/vector specifying time point of interest. (Default: Time-to-event at 50th quantile points)\cr
#' @param newx A numeric/vector specifying biomarker of interest.\cr
#' @param params.x A named vector for biomarker's parameter.\cr
#' @param params.t A named vector for time-to-event's parameter.\cr
#' @param copula A string indicating the type of copula to be used.\cr
#' @param params.copula An integer specifying the copula's parameter.\cr
#' @param params.ph An integer specifying the PH parameter.\cr
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
#' @returns A list of ROC dataframe for each specified time-to-event.
#' @export
timeroc_predict <- function(obj, B = 1, newx, cutoff = 100, t,
                         params.x, params.t, copula,
                         params.copula, params.ph){

  ## preprocessing of arguments
  cov.x <- cov.t <- se.c <- x.val <- dist.x <- dist.t <- ci <- dat <- NULL
  if(inherits(obj, 'fitTROC')){
    args <- preproc(c(as.list(environment()), call = match.call()),
                  extract_from_fitTROC)
  }else if(inherits(obj, 'TimeROC')){
    args <- preproc(c(as.list(environment()), call = match.call()),
                    TimeROC_predict_process)
  }

  list2env(args, environment())

  ### ToDo: Future work for landmark analysis
  # res <- c()
  # for(ttime in t){
  #   subdat <- dat
  #   subdat[which(dat$t > ttime),3] <- 0
  #   fitval <- timeroc_fit(obj, x = subdat$x, t = subdat$t, event = subdat$event)
  #   params.x <- fitval$x$par
  #   cov.x <- fitval$x$cov
  #   if(!missing(params.copula)){
  #     params.t <- fitval$t$par
  #     params.copula <- fitval$copula$par
  #     se.c <- fitval$copula$se
  #     copula <- obj$copula$cop.abbr}
  #   else{
  #     params.t <- fitval$t$par[-length(fitval$t$par)]
  #     params.ph <- fitval$t$par[length(fitval$t$par)]
  #     cov.t <- fitval$t$cov}
  #
  #   pargs <- list(par.x = params.x, par.t = params.t,
  #                 copula = copula,
  #                 par.cop = ifelse((!missing(params.copula)),params.copula,NA),
  #                 par.ph = ifelse((!missing(params.ph)),params.ph,NA),
  #                 cov.x = cov.x, cov.t = cov.t, se.c = se.c,
  #                 t = ttime)
  #   sim.par <- do.call(pars.boot,list(X=x.val, boot.val=B, pargs=pargs))
  #   ret <- do.call(troc.curve,list(pars=sim.par, cutoff=x.val, t=ttime,
  #                                  x.dist = dist.x, t.dist = dist.t,
  #                                  copula = copula, ci = ci))
  #   res <- c(res,ret)
  # }

  pargs <- list(par.x = params.x, par.t = params.t,
                copula = copula,
                par.cop = ifelse((!missing(params.copula)),params.copula,NA),
                par.ph = ifelse((!missing(params.ph)),params.ph,NA),
                cov.x = cov.x, cov.t = cov.t, se.c = se.c,
                t = t)

  sim.par <- do.call(pars.boot,list(X=x.val, boot.val=B, pargs=pargs))
  res <- do.call(troc.curve,list(pars=sim.par, cutoff=x.val, t=t,
                                 x.dist = dist.x, t.dist = dist.t,
                                 copula = copula, ci = ci))
  class(res) <- "predictTROC"
  return(res)
}

# --------------------------------------

#' troc.curve
#'
#' @description helper function to estimate sensitivity and specificity.
#'
#' @param pars A list of vector specifying the parameters for X, T and PH/copula.\cr
#' @param cutoff A vector of integer specifying the biomarker value to be used in the ROC computation.\cr
#' @param t A vector of integer specifying the time to produce the ROC curve.\cr
#' @param x.dist A string to specify the biomarker distribution.\cr
#' @param t.dist A string to specify the time-to-event distribution.\cr
#' @param copula A string to specify the copula distribution.\cr
#' @param ci An integer to specify the confidence interval for the estimated ROC curve.\cr
#' @returns A list of dataframe.
#' @keywords internal
troc.curve <- function(pars, cutoff, t, x.dist, t.dist, copula, ci){
  if(is.null(pars$par.c)) model <- 'ph'
  else model <- 'copula'
  est <- do.call(timeroc,list(pars=pars, cutoff=cutoff, t=t, model=model,
                     x.dist = x.dist, t.dist = t.dist, copula = copula,
                     ci = ci))
  return(est)
}

# --------------------------------------

#' pars.boot
#'
#' @description helper function to generate random parameters used in bootstrap process.
#'
#' @param X A sequence of integer specifying the cutoff to be used.\cr
#' @param boot.val An integer specifying the iteration for the bootstrap process.\cr
#' @param pargs A list of named vector specifying the parameters for X, T and PH/copula.\cr
#' @returns A list of random parameters.
#' @keywords internal
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
pars.boot <- function(X, boot.val, pargs){
  i <- ifelse(is.na(pargs$par.cop),0,1)
  sim <- list()
  if(boot.val <= 1){
    sim$par.x <- matrix(rep(pargs$par.x,2),ncol=length(pargs$par.x), byrow=T)
    colnames(sim$par.x) <- names(pargs$par.x)
    sim$par.t <- matrix(rep(pargs$par.t,2),ncol=length(pargs$par.t), byrow=T)
    colnames(sim$par.t) <- names(pargs$par.t)
    if(i==0){
      sim$par.b <- matrix(rep(pargs$par.ph,2),ncol=length(pargs$par.ph), byrow=T)
      colnames(sim$par.b) <- names(pargs$par.ph)
    } else {
      sim$par.c <- matrix(rep(pargs$par.cop,2),ncol=length(pargs$par.cop), byrow=T)
      colnames(sim$par.c) <- names(pargs$par.cop)
    }

    return(sim)
  } else{
      sim$par.x <- rmvnorm(boot.val, pargs$par.x,pargs$cov.x)
      if(i==0){
        sim$par.t <- rmvnorm(boot.val, c(pargs$par.t,pargs$par.ph), pargs$cov.t)
        sim$par.b <- sim$par.t[,length(ncol(sim$par.t))]
        sim$par.t <- sim$par.t[,-ncol(sim$par.t)]
      } else {
        sim$par.t <- rmvnorm(boot.val, c(pargs$par.t), pargs$cov.t)
        sim$par.c <- rnorm(boot.val, pargs$par.cop, pargs$se.c)
      }
      return(sim)
  }
}

# --------------------------------------

#' timeroc
#'
#' @description helper function that automatically called within troc_curve function.
#'
#' @param model A string of either PH/copula.\cr
#' @param pars A list of X, T and PH/copula parameters.\cr
#' @param t A vector of integer specifying the time to produce the ROC curve.\cr
#' @param cutoff A vector of integer specifying the biomarker value needed for sensitivity & specificity computation.\cr
#' @param x.dist A string specifying the biomarker distribution.\cr
#' @param t.dist A string specifying the time-to-event distribution.\cr
#' @param copula A string specifying the type of copula.\cr
#' @param ci An integer specifying the confidence interval for the estimated ROC curve.\cr
#' @keywords internal
timeroc <- function(model, pars, t, cutoff, x.dist, t.dist, copula, ci){

  xfn <- get.distributions[[x.dist]]$density
  tfn <- get.distributions[[t.dist]]$density
  if(x.dist == "skewnormal") class(xfn) <- "snorm" # there is no sn::psn(..,lower.tail=F)
  if(t.dist == "skewnormal") class(tfn) <- "snorm"
  if(model == 'copula') {
    cfn <- get.copula[[copula]]
    est <- do.call(roc.cop,list(x = cutoff, t = t, xfn = xfn, tfn = tfn,
                                cfn = cfn, pars = pars, ci = ci))
  } else {
    est <- do.call(roc.ph,list(x = cutoff, t = t, xfn = xfn, tfn = tfn,
                               pars = pars, ci = ci))
  }

  return(est)
}

# --------------------------------------

#' roc.ph
#'
#' @description Function that will be automatically called within troc_curve & timeroc to compute sensitivity and specificity of a PH mdoel.
#' @param x An integer vector specifying the biomarker value.\cr
#' @param t An integer vector specifying the time to produce the ROC curve.\cr
#' @param xfn A list of functions used to compute the biomarker distribution.\cr
#' @param tfn A list of functions used to compute the time-to-event distribution.\cr
#' @param pars A list of vector specifying the estimated X, T and PH/copula parameter.\cr
#' @param ci An integer specifying the confidence interval of the ROC curve.\cr
#' @keywords internal
#' @importFrom cubature hcubature
#' @importFrom stats quantile na.omit
roc.ph <- function(x,t, xfn, tfn, pars, ci){
  pobs <- (1-ci)/2
  res <- list()
  res.x <- matrix(NA, nrow = length(x), ncol = 6,
                  dimnames = list(NULL, c('est.sens', 'est.spec', 'low.sens',
                                          'low.spec', 'upp.sens', 'upp.spec')))
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

  rr <- function(args){
    jj <- function(cc) {
      dx <- do.call(xfn[[2]], c(list(x = cc), args[4:length(args)]))
      return(args[2]^(exp(args[3] * cc)) * dx)
    }
    hcubature(jj, lowerLimit = args[1], upperLimit = Inf)$integral
  }

  for(t_val in t){
    pt <- apply(pars$par.t,1,Survt)

    for(ii in seq_along(x)){
      x_val <- x[ii]
      px <- apply(pars$par.x,1,Survx)

      S.t <- apply(cbind(-Inf,pt,pars$par.b,pars$par.x),1,rr)
      s.ct <- apply(cbind(x_val,pt,pars$par.b,pars$par.x),1,rr)

      # sensitivity
      sens <- (px-s.ct) / (1-S.t)
      # sens[which(sens < 0)] <- 0
      # sens[which(sens > 1)] <- 1

      # specificity
      spec <- 1 - (s.ct/S.t)

      res.boot <- na.omit(cbind(sens, spec))
      est <- t(colMeans(res.boot))
      low <- t(apply(res.boot, 2, quantile, probs = pobs))
      upp <- t(apply(res.boot, 2, quantile, probs = 1-pobs))
      res.x[ii,] <- cbind(est, low, upp)
    }
    res[[as.character(t_val)]] <- res.x

  }

  return(res)
}

# --------------------------------------

#' roc.cop
#'
#' @description Function that will be automatically called within troc_curve & timeroc to compute sensitivity and specificity of a copula mdoel.
#' @param x An integer vector specifying the biomarker value.\cr
#' @param t An integer vector specifying the time to produce the ROC curve.\cr
#' @param xfn A list of functions used to compute the biomarker distribution.\cr
#' @param tfn A list of functions used to compute the time-to-event distribution.\cr
#' @param cfn A list of functions used to compute the copula distribution.\cr
#' @param pars A list of vector specifying the estimated X, T and PH/copula parameter.\cr
#' @param ci An integer specifying the confidence interval of the ROC curve.\cr
#' @keywords internal
#' @importFrom stats quantile na.omit
roc.cop <- function(x,t, xfn, tfn, cfn, pars, ci){
  pobs <- (1-ci)/2
  res <- list()
  res.x <- matrix(NA, nrow = length(x), ncol = 6,
                  dimnames = list(NULL, c('est.sens', 'est.spec', 'low.sens',
                                          'low.spec', 'upp.sens', 'upp.spec')))
  res.boot <- matrix(NA, nrow = length(pars$par.b), ncol = 2)
  colnames(res.boot) <- c('sens','spec')

  Ft <- function(pars){
    do.call(tfn[[3]], c(list(t_val),pars))
  }

  Fx <- function(pars){
    do.call(xfn[[3]], c(list(x_val),pars))
  }

  rr <- function(args){
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

      res.boot <- apply(cbind(u,v,pars$par.c),1,rr)
      res.boot <- matrix(unlist(res.boot),ncol=2, nrow=length(pars$par.c),
                         byrow = TRUE)
      est <- na.omit(colMeans(res.boot))
      low <- apply(res.boot, 2, quantile, probs = pobs, na.rm = T)
      upp <- apply(res.boot, 2, quantile, probs = 1-pobs, na.rm = T)
      res.x[ii,] <- cbind(est, low, upp)
    }
    res[[as.character(t_val)]] <- res.x

  }

  return(res)
}
