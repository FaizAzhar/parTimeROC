#' timeroc_auc
#'
#' @description Function to compute the area under the ROC curve
#' @param obj A predictTROC object.\cr
#'
#' @return A dataframe of the area under the ROC curve
#' @export
#' @examples
#' test <- timeroc_obj('normal-weibull-copula', copula = 'clayton90')
#' print(test)
#'
#' set.seed(23456)
#' rr <- rtimeroc(obj = test, censor.rate = 0.1, n=500,
#'                params.t = c(shape=1, scale=5),
#'                params.x = c(mean=5, sd=1),
#'                params.copula=-2)
#'
#' cc <- timeroc_fit(x=rr$x, t=rr$t, event=rr$event, obj = test)
#'
#' jj <- timeroc_predict(cc, t = quantile(rr$t, probs = c(0.25,0.5,0.75)))
#'
#' timeroc_auc(jj)
#' @importFrom DescTools AUC

timeroc_auc <- function(obj){
  if(!inherits(obj, 'predictTROC')) stop("Please provide a predictTROC object")
  res <- data.frame('time' = names(obj),
                    'assoc' = NA,
                    'est.auc' = NA,
                    'low.auc' = NA,
                    'upp.auc' = NA)
  for(i in seq_len(length(obj))){
    res$assoc[i] <- obj[[i]][1,7]
    res$est.auc[i] <- AUC(1-obj[[i]][,2],obj[[i]][,1])
    res$low.auc[i] <- AUC(1-obj[[i]][,4],obj[[i]][,3])
    res$upp.auc[i] <- AUC(1-obj[[i]][,6],obj[[i]][,5])
  }

  return(res)
}
