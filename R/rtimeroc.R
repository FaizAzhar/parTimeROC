#' rtimeroc
#'
#' @description Function to generate bivariate data from PH or copula model.
#'
#' @param obj An initialized 'TimeROC' object.\cr
#' @param n An integer of sample size.\cr
#' @param censor.rate An integer between 0 to 1 that is used for randomized censoring.\cr
#' @param params.x Vector of biomarker parameter.\cr
#' @param params.t Vector of time-to-event parameter.\cr
#' @param params.copula An integer for copula parameter.\cr
#' @param params.ph An integer for association parameter.
#' @export
#' @returns A dataframe with 3 columns (x = biomarker value, t = observable time-to-event, status = censored/not censor (0 or 1))
#' @examples
#' ## Copula model
#' test <- timeroc_obj(dist = 'gompertz-gompertz-copula', copula = "gumbel90")
#' set.seed(23456)
#' rr <- rtimeroc(obj = test, censor.rate = 0, n=500,
#'                params.t = c(shape=3,rate=1),
#'                params.x = c(shape=1,rate=2),
#'                params.copula=-5)
#' plot(t~x, rr)
#'
#' ## PH model
#' test <- timeroc_obj(dist = 'weibull-gompertz-PH')
#' set.seed(23456)
#' rr <- rtimeroc(obj = test, censor.rate = 0, n=100,
#'                params.t = c(shape=2, rate=1),
#'                params.x = c(shape=2, scale=1),
#'                params.ph=0.5)
#' plot(t~x, rr)

rtimeroc <- function(obj, n, censor.rate = 0, params.x, params.t, params.copula, params.ph){

  ## preprocessing of arguments
  copula <- x.dist <- t.dist <- NULL
  args <- preproc(c(as.list(environment()), call = match.call()),
                  extract_from_TimeROC)
  list2env(args, environment())


  if(is.null(as.vector(params.x))) stop("Please provide the parameter for biomarker")
  if(is.null(as.vector(params.t))) stop("Please provide the parameter for time-to-event")
  if(is.na(params.copula) & is.list(copula)) stop("Please provide the parameter for copula")
  if(!is.list(copula) & is.na(params.ph)) stop("Please provide the parameter of PH to be used")
  if(missing(n)) stop("Please provide number of sample to be simulated")


  # generating data
  if(!is.list(copula)){
    rbase <- as.call(c(list(x.dist$density[[1]],n = n), as.list(params.x)))
    x <- eval(rbase)
    S.tx <- stats::runif(n)
    S.t0 <- exp(log(S.tx)/exp(x * params.ph))
    targs <- as.list(params.t)
    if (t.dist$name == "Skew-Normal") qbase <- as.call(c(list(t.dist$density[[4]], p = 1-S.t0), targs))
    else qbase <- as.call(c(list(t.dist$density[[4]], p = S.t0, lower.tail = FALSE), targs))
    t <- eval(qbase)
  } else if(is.list(copula)){
    dat <- as.call(list(copula$density[[1]], N = n, family = copula$family, par = params.copula))
    dd <- eval(dat)
    xargs <- as.list(params.x); targs <- as.list(params.t)
    x <- as.call(c(list(x.dist$density[[4]], p = dd[,1]), xargs))
    t <- as.call(c(list(t.dist$density[[4]], p = dd[,2]), targs))
    x <- eval(x); t <- eval(t)
  }

  status <- rep(1,n)
  tot.c <- floor(censor.rate * n)
  indx <- sample(1:n, tot.c)
  status[indx] <- 0

  return(data.frame(x = x, t = t, event = status))
}
