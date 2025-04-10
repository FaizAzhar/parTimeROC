#
# Store lists of distributions and copulas.
#
# Following flexsurv (C.Jackson), survreg (T.M. Therneau), VineCopula (T. Nagler) framework,
# use closed form hazard h & cumulative hazard H, whenever possible for stable computation.
#
# If do not have a closed formula, h = d/1-F, H = -log((1-F))
#

library(stats)

#' get.distributions
#'
#' Storing list of distributions for Biomarker, X and Time-to-event, T.
#' @returns A list of distributions.
#' @examples
#' get.distributions
#' #> "exponential" "weibull"     "gaussian"    "normal"      "lognormal"
#' #> "gompertz"    "skewnormal"
#' @export
get.distributions <- list(
  'exponential' = list(
    name = "Exponential",
    name.abbr = 'exponential',
    pars = c("rate"),
    init = function(t,...){list(rate = 1/mean(t))},
    density = cbind(rexp, dexp, pexp, qexp),
    cum.hazard = function(t,parms){parms[1] * t},
    hazard = function(t, parms) {parms[1]},
    lower = c(0.01),
    upper = c(1000)
  ),

  'weibull' = list(
    name = "Weibull",
    name.abbr = 'weibull',
    pars = c("shape","scale"),
    init = function(t, x=NULL,event=NULL){
      if (is.null(x) & is.null(event)){
        m <- mean(log(t)); v <- var(log(t))
        tr.scale <- sqrt(v)/1.2; tr.shape <- m + 0.572*tr.scale # following fitdistr package
      }
      else{
        sr.weib <- survival::survreg(survival::Surv(t,event)~x,
                           data=data.frame(x=x,t=t,event=event),
                           dist="weibull")
        tr.shape <- sr.weib$coefficients[1]
        tr.scale <- sr.weib$scale
      }
      list(shape = 1/tr.scale, scale = exp(tr.shape))
    },
    density = cbind(rweibull, dweibull, pweibull, qweibull),
    cum.hazard = function(t, parms){(t/parms[2])^(parms[1])},
    hazard = function(t, parms) {(parms[1]/parms[2]) *
        (t/parms[2])^(parms[1]-1)},
    lower = c(0.1,0.1),
    upper = c(Inf, Inf)
  ),

  'gaussian' = list(
    name = "Gaussian",
    name.abbr = 'normal',
    pars = c("mean","sd"),
    init = function(t, ...){list(mean = mean(t), sd = sd(t))},
    density = cbind(rnorm, dnorm, pnorm, qnorm),
    hazard = function(t, parms) {dnorm(t, mean = parms[1], sd = parms[2])/
        pnorm(t, mean = parms[1], sd = parms[2], lower.tail = FALSE)},
    cum.hazard = function(t, parms) {-pnorm(t, mean = parms[1],
                                            sd = parms[2], lower.tail = FALSE,
                                            log.p = TRUE)},
    lower = c(-Inf, 0.01),
    upper = c(Inf, Inf)
  ),

  'normal' = list(
    name = "Gaussian",
    name.abbr = 'normal',
    pars = c("mean","sd"),
    init = function(t, ...){list(mean = mean(t), sd = sd(t))},
    density = cbind(rnorm, dnorm, pnorm, qnorm),
    hazard = function(t, parms) {dnorm(t, mean = parms[1], sd = parms[2])/
        pnorm(t, mean = parms[1], sd = parms[2], lower.tail = FALSE)},
    cum.hazard = function(t, parms) {-pnorm(t, mean = parms[1],
                                            sd = parms[2], lower.tail = FALSE,
                                            log.p = TRUE)},
    lower = c(-Inf, 0.01),
    upper = c(Inf, Inf)
  ),

  'lognormal' = list(
    name = "Log-Normal",
    name.abbr = 'lognormal',
    pars = c("meanlog","sdlog"),
    init = function(t, x=NULL,event=NULL){
      if (is.null(x) & is.null(event)){
        n <- length(t)
        lt <- log(t)
        mx <- mean(lt)
        sdx <- sqrt((n-1)/n)*sd(lt)
      }
      else{
        sr.lnorm <- survival::survreg(survival::Surv(t,event)~x,
                            data=data.frame(x=x,t=t,event=event),
                            dist="lognormal")
        mx <- sr.lnorm$coefficients[1]
        sdx <- sr.lnorm$scale
      }
      list(meanlog = mx, sdlog = sdx)
    },
    density = cbind(rlnorm, dlnorm, plnorm, qlnorm),
    hazard = function(t,parms){
      logdens <- dlnorm(x = t, meanlog = parms[1],
                        sdlog = parms[2], log = TRUE)
      logsurv <- plnorm(q = t, meanlog = parms[1],
                        sdlog = parms[2], lower.tail = FALSE, log.p = TRUE)
      loghaz <- logdens - logsurv
      exp(loghaz)
    },
    cum.hazard = function(t, parms){
      -plnorm(t, meanlog = parms[1], sdlog = parms[2], lower.tail = FALSE,
              log.p = TRUE)
    },
    lower = c(0.01, 0.01),
    upper = c(Inf, Inf)
  ),

  'gompertz' = list(
    name = "Gompertz",
    name.abbr = 'gompertz',
    pars = c("shape","rate"),
    init = function(t, x=NULL,event=NULL){
      if (is.null(x) & is.null(event)){
        shape <- -0.001
        rate <- 1/mean(t)
      }
      else{
        fsr.gompertz <- flexsurv::flexsurvreg(survival::Surv(t,event)~x,
                                    data=data.frame(x=x,t=t,event=event),
                                    dist="gompertz")
        shape <- fsr.gompertz$res[1,1]
        rate <- fsr.gompertz$res[2,1]
      }
      list(shape = shape, rate = abs(rate))
    },
    density = c(flexsurv::rgompertz, flexsurv::dgompertz, flexsurv::pgompertz,
                flexsurv::qgompertz, flexsurv::hgompertz, flexsurv::Hgompertz),
    hazard = function(t,parms){
      get.distributions[['gompertz']]$density[[5]](t, shape = parms[1], rate = parms[2])
    },
    cum.hazard = function(t, parms){
      get.distributions[['gompertz']]$density[[6]](t, shape = parms[1], rate = parms[2])
    },
    lower = c(-Inf,0.001),
    upper = c(Inf, Inf)
  ),

  'skewnormal' = list(
    name = "Skew-Normal",
    name.abbr = 'skewnormal',
    pars = c("xi","omega",'alpha'),
    init = function(t, ...){
      # Using method-of-moment by Azzalini (because have closed form)
      s <- moments::skewness(t)
      d <- s/sqrt(1+s^2)

      omega <- sqrt(1-(2/pi)*d^2)
      v <- omega^2
      xi <- sqrt(2/pi) * (s/(sqrt(1+s^2)))
      alpha <- ((4-pi)/2)*xi^3/(v^1.5)
      list(xi = xi, omega = omega, alpha = alpha)
    },
    density = c(sn::rsn, sn::dsn, sn::psn, sn::qsn),
    hazard = function(t,parms){
      logdens <- log(sn::dsn(t, xi = parms[1], omega = parms[2], alpha = parms[3]))
      logsurv <- log(1-sn::psn(t, xi = parms[1], omega = parms[2], alpha = parms[3]))
      loghaz <- logdens - logsurv
      exp(loghaz)
    },
    cum.hazard = function(t, parms){
      -log(1-sn::psn(t, xi = parms[1], omega = parms[2], alpha = parms[3]))
    },
    lower = c(-Inf, 0.01, -Inf),
    upper = c(Inf, Inf, Inf)
  ),

  'pch' = list(
    name = "Piecewise Hazard",
    name.abbr = 'pch',
    pars = c("breakpoints","rates"),
    init = function(breakpoints, ...){
      a <- rep(0.1,length(breakpoints)-1)
      list(breakpoints = breakpoints, rates = a)
    },
    density = c(rpch = function(n, breakpoints, rates){
                  p <- runif(n)
                  get.distributions$pch$density$qpch(p, breakpoints, rates)
                },
                dpch = function(x, breakpoints, rates, lower.tail = TRUE, log.p = FALSE, ...){
                  res <- get.distributions$pch$hazard(x, breakpoints, rates,...) *
                    exp(-get.distributions$pch$cum.hazard(x, breakpoints, rates,...))
                  res
                },
                ppch = function(q, breakpoints, rates, lower.tail = TRUE, log.p = FALSE, ...){
                  H_t <- get.distributions$pch$cum.hazard(q, breakpoints, rates,...)
                  res <- 1 - exp(-H_t)
                  if (log.p & lower.tail){
                    res <- log(res)
                  } else if(log.p & !lower.tail){
                    res <- log(1-res)
                  } else if(!log.p & !lower.tail){
                    res <- 1-res
                  }
                  res
                },
                qpch = function(p, breakpoints, rates, lower.tail = TRUE, log.p = FALSE, ...){
                  if(length(rates) != length(breakpoints) - 1) stop("Length rated is not equal to total intervals")
                  if (log.p) p <- exp(p)
                  if(any(p >= 1)) stop("Quantile must be < 1")
                  if(any(p <= 0)) stop("Quantile must be > 0")
                  if(!lower.tail) p <- 1-p
                  res <- numeric(0)
                  H_t <- c(0,cumsum(diff(breakpoints) * rates))

                  for(i in seq_along(p)){
                    Hcheck <- -log(1-p[i])
                    pos <- findInterval(Hcheck, vec = H_t)
                    if(breakpoints[pos+1] == max(breakpoints) & !is.infinite(breakpoints[pos+1])) upper <- breakpoints[pos+1] - 0.01
                    else if(is.numeric(breakpoints[pos+1])) upper <- breakpoints[pos+1]
                    res[i] <- uniroot(f=function(q){
                      get.distributions$pch$density$ppch(q,breakpoints,rates) - p[i]
                    }, interval = c(breakpoints[pos],upper))$root
                  }
                  res
                }),
    hazard = function(t, breakpoints, rates,...){
      setup <- setup.pch(t, breakpoints)
      h <- setup$mat_event %*% rates
      h
    },
    cum.hazard = function(t, breakpoints, rates,...){
      setup <- setup.pch(t, breakpoints)
      H <- setup$mat_exposure %*% rates
      H
    }
  ),

  'llogis' = list(
    name = "Log-logistic",
    name.abbr = 'llogis',
    pars = c("shape","scale"),
    init = function(t, ...){
      list(shape = 1, scale = 1)
    },
    density = c(flexsurv::rllogis, flexsurv::dllogis, flexsurv::pllogis, flexsurv::qllogis),
    hazard = function(t,parms){
      logdens <- log(flexsurv::dllogis(t, shape = parms[1], scale = parms[2]))
      logsurv <- log(1-flexsurv::pllogis(t, shape = parms[1], scale = parms[2]))
      loghaz <- logdens - logsurv
      exp(loghaz)
    },
    cum.hazard = function(t, parms){
      -log(1-flexsurv::pllogis(t, shape = parms[1], scale = parms[2]))
    },
    lower = c(0, 0),
    upper = c(Inf, Inf)
  )
)

#' get.copula
#'
#' Storing list of copula.
#' @returns A list of copula.
#' #' @examples
#' get.copula
#' #> "gaussian"  "clayton90" "gumbel90"  "gumbel"    "joe90"
#' @export
get.copula <- list(
  'gaussian' = list(
    name = "Gaussian",
    cop.abbr = "gaussian",
    family = 1,
    init = function(u,v,...){
      tau <- VineCopula::BiCopEst(u,v, family = 1)$emptau
      VineCopula::BiCopTau2Par(family = 1, tau)
      },
    density = cbind(VineCopula::BiCopSim, VineCopula::BiCopPDF,
                    VineCopula::BiCopCDF, VineCopula::BiCopHfunc1,
                    VineCopula::BiCopHfunc2),
    lower = -0.9999999,
    upper = 0.9999999
  ),

  'clayton90' = list(
    name = "90 Degrees Rotated Clayton",
    cop.abbr = "clayton90",
    family = 23,
    init = function(u,v,...){
      tau <- VineCopula::BiCopEst(u,v, family = 23)$emptau
      VineCopula::BiCopTau2Par(family = 23, tau)
      },
    density = cbind(VineCopula::BiCopSim, VineCopula::BiCopPDF,
                    VineCopula::BiCopCDF, VineCopula::BiCopHfunc1,
                    VineCopula::BiCopHfunc2),
    lower = -27.99,
    upper = -0.01
  ),

  'gumbel90' = list(
    name = "90 Degrees Rotated Gumbel",
    cop.abbr = "gumbel90",
    family = 24,
    init = function(u,v,...){
      tau <- VineCopula::BiCopEst(u,v, family = 24)$emptau
      VineCopula::BiCopTau2Par(family = 24, tau)
      },
    density = cbind(VineCopula::BiCopSim, VineCopula::BiCopPDF,
                    VineCopula::BiCopCDF, VineCopula::BiCopHfunc1,
                    VineCopula::BiCopHfunc2),
    lower = -16.99,
    upper = -1.01
  ),

  'gumbel' = list(
    name = "Gumbel",
    cop.abbr = "gumbel",
    family = 4,
    init = function(u,v,...){
      tau <- VineCopula::BiCopEst(u,v, family = 4)$emptau
      VineCopula::BiCopTau2Par(family = 4, tau)
    },
    density = cbind(VineCopula::BiCopSim, VineCopula::BiCopPDF,
                    VineCopula::BiCopCDF, VineCopula::BiCopHfunc1,
                    VineCopula::BiCopHfunc2),
    lower = 16.99,
    upper = 1.01
  ),

  'joe90' = list(
    name = "90 Degrees Rotated Joe",
    cop.abbr = "joe90",
    family = 26,
    init = function(u,v,...){
      tau <- VineCopula::BiCopEst(u,v, family = 26)$emptau
      VineCopula::BiCopTau2Par(family = 26, tau)
      },
    density = cbind(VineCopula::BiCopSim, VineCopula::BiCopPDF,
                    VineCopula::BiCopCDF, VineCopula::BiCopHfunc1,
                    VineCopula::BiCopHfunc2),
    lower = -29.99,
    upper = -1.01
  )
)
