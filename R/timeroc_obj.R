#' timeroc_obj
#'
#' @description Function to initialized time-dependent ROC object.
#'
#' @param dist A string emphasizing the distribution assumption for biomarker-time-model.\cr
#' @param params.x params.x Vector of biomarker parameter.\cr
#' @param params.x Vector of biomarker parameter.\cr
#' @param params.t Vector of time-to-event parameter.\cr
#' @param copula A string emphasizing on the copula to be used.\cr
#' @param params.copula An integer for copula parameter.\cr
#' @param params.ph An integer for association parameter.
#' @export
#' @import flexsurv
#' @import survival
#' @import VineCopula
#' @import sn
#' @returns A 'TimeROC' object.
#' @examples
#' ## Copula model
#' test <- timeroc_obj(dist = 'gompertz-gompertz-copula', copula = "gumbel90")
#'
#' ## PH model
#' test <- timeroc_obj(dist = 'weibull-gompertz-PH')

timeroc_obj <- function(dist, params.x = NA, params.t = NA, copula = NA,
                        params.copula = NA, params.ph = NA){

  x.dist = NA; t.dist = NA
  if (missing(dist)){stop("Please provide the distributions and models to be used.
                          EG: normal-weibull-PH / normal-weibull-copula")}
  else {
    dlist <- unlist(strsplit(tolower(dist), split = "-"))

    # check list distribution & copula
    if(dlist[1] %in% names(get.distributions)){x.dist <- get.distributions[[dlist[1]]]}
    else{stop("Biomarker distribution not in the list")}
    if(dlist[2] %in% names(get.distributions)){t.dist <- get.distributions[[dlist[2]]]}
    else{stop("Time-to-event distribution not in the list")}
    if(dlist[3] == 'copula' & is.na(copula)) stop("Please provide a listed copula")
    else if(dlist[3] == 'copula' & !is.na(copula)){
      copula <- get.copula[[copula]]
      if(is.null(copula)) stop("Copula not exist in list")
    }
  }

  obj <- list(
    dist = dist,
    x.dist = x.dist,
    t.dist = t.dist,
    copula = copula,
    params.x = params.x,
    params.t = params.t,
    params.copula = params.copula,
    params.ph = params.ph
  )

  class(obj) <- "TimeROC"
  return(obj)
}
