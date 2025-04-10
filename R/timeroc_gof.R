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
#' @importFrom GofCens KScens

timeroc_gof <- function(obj){
  ## preprocessing of arguments
  params.t <- params.x <- copula <- params.ph <- params.copula <- NULL
  iscopula <- x.dist <- t.dist <- breakpoints<- NULL
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

    # # bootstrap Kolmogorov-Smirnov test
    # boot_obj <- timeroc_obj(paste0(x.dist$name.abbr,"-",t.dist$name.abbr, "-PH"))
    # NN <- nrow(df)
    # ks_boot <- c()
    # for(i in 1:1000){
    #   idx <- sample(1:NN, NN, replace = TRUE)
    #   boot_sam <- rtimeroc(boot_obj, n = NN,
    #                        params.x = params.x,
    #                        params.t = params.t,
    #                        params.ph = params.ph)
    #   boot_fit <- timeroc_fit(boot_obj, x = boot_sam$x, t = boot_sam$t, event = boot_sam$event)
    #   boot_coxsnell <- t.dist$cum.hazard(boot_sam$t, boot_fit$t$par) * exp(boot_sam$x * boot_fit$t$par[length(boot_fit$t$par)])
    #   # EDF-test based on Cockeran et al., (2019)
    #   order_coxsnell <- data.frame('order' = 1:NN, 'coxsnell' = sort(boot_coxsnell))
    #   ks_pos <- max(order_coxsnell$order/NN-(1-exp(-order_coxsnell$coxsnell)))
    #   ks_neg <- max((1-exp(-order_coxsnell$coxsnell))-(order_coxsnell$order - 1)/NN)
    #   ks_boot <- append(ks_boot,max(ks_pos,ks_neg))
    # }
    #
    # ks_dist <- ecdf(ks_boot)

    # Cox-Snell Residual
    df <- df[order(df$x),]
    df <- df[order(df$t),]
    if(t.dist$name.abbr != "pch"){
      df$coxsnell <- t.dist$cum.hazard(df$t, params.t) * exp(df$x * params.ph)
    } else{
        df$coxsnell <- as.vector(t.dist$cum.hazard(df$t, breakpoints, rates=params.t) * exp(df$x * params.ph))
    }
    # df$coxsnell[which(df$event == 0)] <- df$coxsnell[which(df$event == 0)] + log(2)
    df$mresid <- df$event - df$coxsnell
    df$sgn <- 1; df$sgn[which(df$mresid < 0)] <- -1
    df$devresid <- df$sgn*sqrt(-2*(df$mresid + df$event * log(df$event - df$mresid)))
    fit.err <- coxph(Surv(coxsnell,event)~1, data = df, method = 'breslow')
    df$Hcoxsnell <- predict(fit.err, type = 'expected')

    # order_coxsnell <- data.frame('order' = 1:NN, 'coxsnell' = sort(df$coxsnell))
    # ks_pos <- max(order_coxsnell$order/NN-(1-exp(-order_coxsnell$coxsnell)))
    # ks_neg <- max((1-exp(-order_coxsnell$coxsnell))-(order_coxsnell$order - 1)/NN)
    # ks_n <- max(ks_pos,ks_neg)
    #
    # # Bolshev correction
    # ks <- (6*NN * ks_n + 1) / (6*sqrt(NN))
    # ks_t <- list(statistic = ks, p.value = 1-ks_dist(ks))

    # # TO-DO (Breslow Estimator & LOESS)
    # fit <- coxph(formula = Surv(t, event) ~ 1, data = df, method = "breslow")
    # surv_table <- basehaz(fit)
    # # naest1 <- cumsum(fit$n.event/fit$n.risk)
    # # surv_table <- data.frame("time" = fit$time, "hazard" = naest1)
    #
    # df <- merge(df, surv_table[, c("time", "hazard")], by.x = "t", by.y = "time", all.x = TRUE)
    # df <- df[order(df$t), ]  # Reordering if necessary
    # colnames(df)[colnames(df) == "hazard"] <- "H_NA"

    # theoretical.q <- ecdf(df$coxsnell)
    # theo.q <- theoretical.q(df$coxsnell)
    # theo.expo <- pexp(df$coxsnell)
    # ks_t <- ks.test(df$coxsnell, 'pexp')
    # ks_t <- ks.test(df$coxsnell, df$Hcoxsnell)
    ks_t <- GofCens::KScens(Surv(coxsnell,event)~1, data = df, distr = "exponential")$Test
    DescTools::PlotQQ(df$coxsnell[which(df$event == 1)], function(p){qexp(p,rate=1)},
      main = paste0("Time-to-event (K-S : p = ",round(ks_t[2],4),")"),
      xlab = "Theoretical Exponential",
      ylab = "Cox-Snell Residuals")

    # df <- df[order(df$coxsnell),]
    # loess_PH <- loess(H_NA~coxsnell, data = df, span = 0.5)
    # smooth_residual <- predict(loess_PH)
    #
    # plot(df$coxsnell, df$H_NA, main = paste0('Time-to-event ',
    #                                            '(K-S : p = ',round(ks_t$p.value,4),
    #                                            ')'),
    #      xlab = "Cox-Snell Residual",
    #      ylab = "Cum.Haz",
    #      xlim = c(0,max(df$coxsnell)),
    #      ylim = c(0,max(df$H_NA)))
    # lines(smooth_residual, x = df$coxsnell, col = 'red')
    # abline(a=0,b=1, lty='dashed', col= 'blue')

    layout(matrix(c(1), nrow = 1, ncol = 1))
    return(list(ks_x = ks_x, ks_t = ks_t, df = df))
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
