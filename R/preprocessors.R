#
# Store functions that will extract information
# from the initialized 'TimeROC' & 'fitTROC' object.
#

preproc <- function(args, ...){
  # what has to be done?
  tasks <- list(...)

  # perform all tasks sequentially
  for (i in seq_along(tasks)) {
    stopifnot(is.function(tasks[[i]]))
    args <- tasks[[i]](args)
  }

  # return preprocessed arguments
  args
}

extract_from_TimeROC <- function(args){
  if (!inherits(args$obj, 'TimeROC')) stop("Please supply a TimeROC object")

  # set to NA if missing from input
  if (is.symbol(args$params.x)) args$params.x <- NA
  if (is.symbol(args$params.t)) args$params.t <- NA
  if (is.symbol(args$params.copula)) args$params.copula <- NA
  if (is.symbol(args$params.ph)) args$params.ph <- NA

  # store info from TimeROC object into environment
  args$dist <- args$obj$dist
  args$x.dist <- args$obj$x.dist
  args$t.dist <- args$obj$t.dist
  args$copula <- args$obj$copula
  if (!all(is.na(args$obj$params.x)))args$params.x <- args$obj$params.x
  if (!all(is.na(args$obj$params.t)))args$params.t <- args$obj$params.t
  if (!is.na(args$obj$params.copula))args$params.copula <- args$obj$params.copula
  if (!is.na(args$obj$params.ph))args$params.ph <- args$obj$params.ph

  args
}

TimeROC_predict_process <- function(args){

  # set to NA if missing from input
  if (is.symbol(args$params.x)) args$params.x <- NA
  if (is.symbol(args$params.t)) args$params.t <- NA
  if (is.symbol(args$params.copula)) args$params.copula <- NA
  if (is.symbol(args$params.ph)) args$params.ph <- NA
  if (is.symbol(args$seed)) args$seed <- NA

  # store info from TimeROC object into environment
  if(!is.list(args$obj$copula)) {
    args$name.copula <- NA
    args$iscopula <- 0
  } else {
    args$name.copula <- args$obj$copula$cop.abbr
    args$copula <- args$obj$copula
    args$iscopula <- 1
  }
  if (!all(is.na(args$obj$params.x)))args$params.x <- args$obj$params.x
  if (!all(is.na(args$obj$params.t)))args$params.t <- args$obj$params.t
  if (!is.na(args$obj$params.copula))args$params.copula <- args$obj$params.copula
  if (!is.na(args$obj$params.ph))args$params.ph <- args$obj$params.ph
  if (!is.na(args$seed))args$seed <- args$seed

  args$name.dist.x <- args$obj$x.dist$name.abbr
  args$name.dist.t <- args$obj$t.dist$name.abbr
  args$x.dist <- args$obj$x.dist
  args$t.dist <- args$obj$t.dist
  args$cum.haz <- args$obj$t.dist$cum.hazard

  if(!is.na(args$seed)) {set.seed(args$seed)}
  rr <- rtimeroc(args$obj, n = 500, censor.rate = 0,
                   params.x = args$params.x,
                   params.t = args$params.t,
                   params.copula = args$params.copula,
                   params.ph = args$params.ph)

  args$xval <- if(is.symbol(args$newx)) {
    seq(min(rr[,1]),max(rr[,1]),length.out = args$cutoff)
  } else{args$newx}

  if(is.symbol(args$t)) args$t <- stats::quantile(rr[,2],probs = 0.5)

  args$se.x <- rep(1,length(args$params.x))
  names(args$se.x) <- names(args$params.x)
  args$se.t <- rep(1,length(args$params.t))
  names(args$se.t) <- names(args$params.t)
  if(is.null(args$params.copula)) args$se.ph <- 1
  else args$se.ph <- NA
  if(!is.null(args$params.copula)) args$se.c <- 1
  else args$se.c <- NA

  ## ToDo: Future work for landmark analysis
  args$dat <- rr

  args
}

extract_from_fitTROC <- function(args){
  if (!inherits(args$obj, 'fitTROC')) stop("Please supply a fitTROC object")

  # set TimeROC object and store fitted parameters
  if(is.null(args$obj$name.copula)) {
    cname <- NA
    args$iscopula <- 0
  } else {
    cname <- args$obj$name.copula
    args$iscopula <- 1
  }
  j <- !(names(args$obj$t$par) %in% 'beta')
  out <- timeroc_obj(args$obj$name, params.x = args$obj$x$par,
                   params.t = args$obj$t$par[j],
                   params.copula = args$obj$copula$par['theta'],
                   params.ph = args$obj$t$par['beta'],
                   copula = cname)

  # store info from fitTROC object into environment
  args$cum.haz <- out$t.dist$cum.hazard
  args$params.x <- out$params.x
  args$params.t <- out$params.t
  args$name.copula <- ifelse(is.null(args$obj$name.copula),NA,out$copula$cop.abbr)
  args$name.dist.x <- out$x.dist$name.abbr
  args$name.dist.t <- out$t.dist$name.abbr
  args$copula <- out$copula
  args$x.dist <- out$x.dist
  args$t.dist <- out$t.dist
  args$params.copula <- out$params.copula
  args$params.ph <- out$params.ph
  args$se.x <- args$obj$x$se
  args$se.t <- args$obj$t$se[j]
  args$xval <- if(is.symbol(args$newx)) {seq(min(args$obj$dat[,1]),
                                                max(args$obj$dat[,1]),
                                                length.out = args$cutoff)} else{args$newx}

  if(is.symbol(args$t)) args$t <- stats::quantile(args$obj$dat[,2],
                                              probs = 0.5)

  if(is.null(args$params.copula)) args$se.ph <- args$obj$t$se['beta']
  else args$se.ph <- NA
  if(!is.null(args$params.copula)) args$se.c <- args$obj$copula$se
  else args$se.c <- NA
  if(!is.null(args$obj$t$breakpoints)) args$breakpoints <- args$obj$t$breakpoints
  else args$breakpoints <- NA

  ## ToDo: Future work for landmark analysis
  args$dat <- args$obj$dat
  # if(!is.null(args$params.copula)){
  #   args$obj <- timeroc_obj(dist = paste(args$dist.x,args$dist.t,'copula', sep = '-'),
  #                           copula = cname)
  # }else{
  #   args$obj <- timeroc_obj(dist = paste(args$dist.x,args$dist.t,'PH', sep = '-'))
  # }

  args
}
