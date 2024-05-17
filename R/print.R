#' @export
print.TimeROC <- function(x, ...){
  cat("Model Assumptions: ")
  if(is.list(x$copula)){
    cat(x$copula$name, "Copula")
  } else {cat("Proportional Hazard (PH)")}

  cat('\nX                :',x$x.dist$name)
  cat('\nTime-to-Event    :',x$t.dist$name,"\n")
}

#' @export
print.fitTROC <- function(x, ...){
  if(is.null(x$copula)){
    tparname <- !(names(x$t$par) %in% c("beta"))
    b.dat <- data.frame(est = round(x$t$par["beta"],4),
                        low = round(x$t$lbound["beta"],4),
                        upper = round(x$t$ubound["beta"],4),
                        se = round(x$t$se["beta"],4))
  } else tparname <- names(x$t$par)

  x.dat <- data.frame(est = round(x$x$par,4),
                      low = round(x$x$lbound,4),
                      upper = round(x$x$ubound,4),
                      se = round(x$x$se,4))

  t.dat <- data.frame(est = round(x$t$par[tparname],4),
                      low = round(x$t$lbound[tparname],4),
                      upper = round(x$t$ubound[tparname],4),
                      se = round(x$t$se[tparname],4))

  if(!is.null(x$copula)){
    c.dat <- data.frame(est = round(x$copula$par,4),
                        low = round(x$copula$lbound,4),
                        upper = round(x$copula$ubound,4),
                        se = round(x$copula$se,4))
  }

  cat('Model: ', x$name,'\n')
  cat('------\n')
  cat('X', paste0('(',x$conf*100,'% CI)'),':\nAIC = ', x$x$aic,'\n')
  print(x.dat)

  cat('------\n')
  cat('Time-to-Event', paste0('(',x$conf*100,'% CI)'),':\nAIC = ', x$t$aic,'\n')
  print(t.dat)
  cat('------\n')
  if(!is.null(x$copula)){cat('Copula', paste0('(',x$conf*100,'% CI)'),':\nAIC = ', x$copula$aic,'\n');
    print(c.dat)}
  else{cat('PH', paste0('(',x$conf*100,'% CI)'),':\n'); print(b.dat)}
}
