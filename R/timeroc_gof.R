#' timeroc_gof
#'
#' @description Function to compute goodness-of-fit for the Proportional Hazard (PH) or copula model. For PH model, the Cox-Snell residuals are computed and compared with Exponential(rate=1). For copula model, the Rosenblatt transformation is applied before performing independence testing. Kolmogorov-Smirnov test is performed to check the goodness-of-fit of the biomarker and time-to-event.
#'
#' @param obj A 'fitTROC' object returned from fitting procedure.
#'
#' @return A list of test statistics and p-values. Automatically plot residuals for biomarker and time-to-event.
#' @export
#'
#' @examples
#' # Copula model
#' rt <- timeroc_obj("normal-weibull-copula",copula="clayton90")
#' set.seed(1)
#' rr <- rtimeroc(rt, n=300, censor.rate = 0,
#'                params.x = c(mean=5, sd=1),
#'                params.t = c(shape=1, scale=5),
#'                params.copula = -2.5)
#' plot(t~x, data=rr)
#' test <- timeroc_obj("normal-weibull-copula",copula="gumbel90")
#' jj <- timeroc_fit(test, rr$x, rr$t, rr$event)
#'
#' cc <- timeroc_gof(jj)
#'
#' test <- timeroc_obj("normal-weibull-copula",copula="clayton90")
#' jj <- timeroc_fit(test, rr$x, rr$t, rr$event)
#'
#' cc <- timeroc_gof(jj)
#'
#' # PH model
#' rt <- timeroc_obj("normal-weibull-PH")
#' set.seed(1)
#' rr <- rtimeroc(rt, n=300, censor.rate = 0,
#'                params.x = c(mean=5, sd=1),
#'                params.t = c(shape=1, scale=5),
#'                params.ph = 1.2)
#' plot(t~x, data=rr)
#' test <- timeroc_obj("lognormal-lognormal-PH")
#' jj <- timeroc_fit(test, rr$x, rr$t, rr$event)
#'
#' cc <- timeroc_gof(jj)
#'
#' test <- timeroc_obj("normal-weibull-PH")
#' jj <- timeroc_fit(test, rr$x, rr$t, rr$event)
#'
#' cc <- timeroc_gof(jj)
#'
#' @importFrom stats ppoints ks.test qexp
#' @importFrom graphics layout

timeroc_gof <- function(obj){
  ## preprocessing of arguments
  params.t <- params.x <- copula <- params.ph <- params.copula <- NULL
  iscopula <- x.dist <- t.dist <- NULL
  df <- obj$dat

  args <- preproc(c(as.list(environment()), call = match.call()),
                  extract_from_fitTROC)
  list2env(args, environment())

  DescTools::DescToolsOptions(stamp=NULL)

  if(!iscopula){
    layout(matrix(c(1, 2), nrow = 1, ncol = 2))

    # QQ plot X
    pp <- ppoints(nrow(df))
    theoretical.q <- as.call(c(list(x.dist$density[[4]], p = pp), as.list(params.x)))
    theo.q <- eval(theoretical.q)
    ks_x <- ks.test(df$x,theo.q)
    DescTools::PlotQQ(df$x, function(p){
      qq <- as.call(c(list(x.dist$density[[4]],p=p),params.x))
      eval(qq)},
      main = paste0("Biomarker (K-S : p = ",round(ks_x$p.value,4),")"),
      xlab = paste0("Theoretical ",x.dist$name))

    # Cox-Snell Residual
    df <- df[order(df$t),]
    df$coxsnell <- t.dist$cum.hazard(df$t, params.t) * exp(df$x * params.ph)
    # H_nelson <- nelson_aalen(df)
    # df$nelson <- H_nelson$Hj[H_nelson$tj == df$t]
    theoretical.q <- call("qexp", p = pp, rate = 1)
    theo.q <- eval(theoretical.q)
    ks_t <- ks.test(df$coxsnell, theo.q)
    DescTools::PlotQQ(df$coxsnell, function(p){qexp(p,rate=1)},
      main = paste0("Time-to-event (K-S : p = ",round(ks_t$p.value,4),")"),
      xlab = "Theoretical Exponential",
      ylab = "Cox-Snell Residuals")

    # plot(nelson~coxsnell, data = df, main = paste0('Time-to-event ',
    #                                            '(K-S : p = ',round(ks_t$p.value,4),
    #                                            ')'),
    #      xlab = "Cox-Snell Residual",
    #      ylab = "Cum.Haz (Nelson-Aalen)")
    # abline(a=0,b=1, lty='dashed', col= 'blue')

    layout(matrix(c(1), nrow = 1, ncol = 1))
    return(list(ks_x = ks_x, ks_t = ks_t))
  }
  else{
    # layout(matrix(c(1, 0, 1,  3, 2, 3, 2, 0), nrow = 2, ncol = 4))
    layout(matrix(c(1, 3, 1,3,2,  4,2,4), nrow = 2, ncol = 4))

    # QQ plot X
    pp <- ppoints(nrow(df))

    theoretical.q <- as.call(c(list(x.dist$density[[4]], p = pp), as.list(params.x)))
    theo.q <- eval(theoretical.q)

    ks_x <- ks.test(df$x,theo.q)
    DescTools::PlotQQ(df$x, function(p){
      qq <- as.call(c(list(x.dist$density[[4]],p=p),params.x))
      eval(qq)},
      main = paste0("Biomarker (K-S : p = ",round(ks_x$p.value,4),")"),
      xlab = paste0("Theoretical ",x.dist$name))

    # QQ plot T
    pp <- ppoints(nrow(df))
    theoretical.q <- as.call(c(list(t.dist$density[[4]], p = pp), as.list(params.t)))
    theo.q <- eval(theoretical.q)
    ks_t <- ks.test(df$t,theo.q)
    DescTools::PlotQQ(df$t, function(p){
      qq <- as.call(c(list(t.dist$density[[4]],p=p),params.t))
      eval(qq)},
      main = paste0("Time-to-event (K-S : p = ",round(ks_t$p.value,4),")"),
      xlab = paste0("Theoretical ",t.dist$name))

    # QQ plot Copula (Using Rosenblatt transformations)
    Ft <- as.call(c(list(t.dist$density[[3]],df$t),params.t))
    Fx <- as.call(c(list(x.dist$density[[3]],df$x),params.x))
    u <- eval(Fx); v <- eval(Ft)
    # Ft <- stats::ecdf(df$t)
    # Fx <- stats::ecdf(df$x)
    # u <- Fx(df$x); v <- Ft(df$t)
    cond_uv <- VineCopula::BiCopHfunc(u,v, family = copula$family, par = params.copula)
    # kendall_test <- VineCopula::BiCopGofTest(u,v, family = cdist$family,
    #                                        par = params.copula,
    #                                        method = "kendall", B = 200)


    # v|u plot
    v_cond <- cond_uv$hfunc1
    indep_u <- VineCopula::BiCopIndTest(u,v_cond)
    plot(x = u, y = v_cond, main = paste0("Pairwise Rosenblatt (Ind.Test : p = ",
                                    round(indep_u$p.value,4),")"),
         xlab = paste0("Empirical u"),
         ylab = "Empirical v|u")
    # u|v plot
    u_cond <- cond_uv$hfunc2
    indep_v <- VineCopula::BiCopIndTest(v,u_cond)
    # yax <- stats::qnorm(u_cond)^2 + stats::qnorm(v)^2
    # DescTools::PlotQQ(yax, function(p){stats::qchisq(p,df=2)},
    #                   main = paste0("Pairwise Rosenblatt (Ind.Test : p = ",
    #                                 round(indep_test$p.value,4),")"),
    #                   xlab = paste0("Theoretical Chi-Square (df = 2)"))
    plot(x = v, y = u_cond, main = paste0("Pairwise Rosenblatt (Ind.Test : p = ",
                                          round(indep_v$p.value,4),")"),
         xlab = paste0("Empirical v"),
         ylab = "Empirical u|v")

    layout(matrix(c(1), nrow = 1, ncol = 1))

    return(list(ks_x = ks_x, ks_t = ks_t, ind_u = indep_u, ind_v = indep_v))
  }

}
