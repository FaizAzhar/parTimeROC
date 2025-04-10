ipcw <- function(df, obj){
  df <- df[order(df$t, df$x, decreasing = FALSE),]
  df$id <- 1:nrow(df)
  df$censored <- 1 - df$event
  # times <- sort(df$t, decreasing = FALSE)
  # test.long <- survSplit(df, cut=times, end="t",
  #                        start="Tstart", event="event")
  # test.long <- test.long[order(test.long$id,test.long$t),]
  # test.long.cens <- survSplit(df, cut = times, end="t",
  #                             start="Tstart", event="censored")
  # test.long.cens <- test.long.cens[order(test.long.cens$id,
  #                                        test.long.cens$t),]
  # test.long$censored <- test.long.cens$censored
  # C0 <- coxph(Surv(Tstart, t, censored)~1, data = test.long)
  # CZ <- coxph(Surv(Tstart, t, censored)~x, data = test.long)
  # C0fit <- summary(survfit(C0), times = test.long$Tstart)
  # test.long$K0ti <- C0fit$surv
  # test.long$KZti <- NULL
  # for(i in 1:nrow(test.long)){
  #   datai <- test.long[i,]
  #   sfiCZ <- survfit(CZ, newdata = datai)
  #   ssfiCZ <- summary(sfiCZ, times = datai$Tstart)
  #   test.long$KZti[i] <- ssfiCZ$surv
  # }
  # test.long$WUnStab <- 1/test.long$KZti
  # test.long$WStab <- test.long$K0ti/test.long$KZti
  # res <- test.long[order(test.long$id, test.long$t, decreasing = TRUE), ]
  # res <- res[!duplicated(res$id), ]

  # tryCatch({
  #   fitC0 <- timeroc_fit(obj, x = df.cens$x, t = df.cens$t, event = df.cens$censored)
  #   }, error = function(e){})

  # cox.mod <- coxph(Surv(t,censored)~x,data=df)
  # cox.mod0 <- coxph(Surv(t,censored)~0,data=df)
  # Sti <- predict(cox.mod, newdata=df.cens, type = 'survival', se.fit=TRUE)
  # Sti0 <- survfit(cox.mod0, newdata=df.cens)
  # df$KZti <- Sti$fit
  # surv_tab <- data.frame(surv=Sti0$surv, time=Sti0$time)
  #
  # df$K0ti <- NULL
  # for (i in 1:nrow(df)){
  #   df$K0ti[i] <- surv_tab[which(surv_tab$time == df$t[i]), 'surv']
  # }


  res.t <- fit.parmsurvival(df$t, df$censored, obj$t.dist)

  df$K0ti <- exp(-obj$t.dist$cum.hazard(df$t, res.t$par))

  df$W <- 1/(df$K0ti)
  df
}

fit.parmsurvival <- function(t, censored, t.dist){
  Call2 <- match.call(expand.dots = TRUE)
  res.t <- list()

  Call2[[1L]] <- quote(stats::optim)
  ll.t <- function(parms, ...){
    likeli <- (censored * (log(t.dist$hazard(t,parms))) -
                         t.dist$cum.hazard(t, parms))
    -sum(likeli)
  }
  Call2$fn <- ll.t
  init.param.t <- t.dist$init(t)
  Call2$par <- unlist(init.param.t)
  Call2$hessian <- TRUE
  Call2$lower <- c(t.dist$lower, -Inf)
  Call2$upper <- c(t.dist$upper, Inf)
  Call2$method <- "L-BFGS-B"
  res.t <- eval.parent(Call2)
  names(res.t$par) <- c(t.dist$pars)
  return(res.t)
}
