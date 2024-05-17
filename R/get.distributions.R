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
      list(shape = shape, rate = rate)
    },
    density = c(flexsurv::rgompertz, flexsurv::dgompertz, flexsurv::pgompertz,
                flexsurv::qgompertz, flexsurv::hgompertz, flexsurv::Hgompertz),
    hazard = function(t,parms){
      get.distributions[['gompertz']]$density[[5]](t, shape = parms[1], rate = parms[2])
    },
    cum.hazard = function(t, parms){
      get.distributions[['gompertz']]$density[[6]](t, shape = parms[1], rate = parms[2])
    },
    lower = c(-Inf,0.01),
    upper = c(Inf, Inf)
  ),

  'skewnormal' = list(
    name = "Skew-Normal",
    name.abbr = 'skewnormal',
    pars = c("xi","omega",'alpha'),
    init = function(t, ...){
      # Using method-of-moment by Azzalini (because have closed form)
      m <- mean(t); v <- sd(t)
      s <- mean((t-m)^3)/(mean((t-m)^2)^(3/2))
      d <- (pi/2)*(abs(s)^(2/3))/(abs(s)^(2/3) + ((4-pi)/2)^(2/3))
      alpha <- d/sqrt(1-d^2)
      omega <- v/sqrt(1-(2*d^2)/pi)
      xi <- m - omega*d*sqrt(2/pi)
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
  )
)

#' get.copula
#'
#' Storing list of copula.
#' @returns A list of copula.
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
                    VineCopula::BiCopCDF, VineCopula::BiCopHfunc1),
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
                    VineCopula::BiCopCDF, VineCopula::BiCopHfunc1),
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
                    VineCopula::BiCopCDF, VineCopula::BiCopHfunc1),
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
                    VineCopula::BiCopCDF, VineCopula::BiCopHfunc1),
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
                    VineCopula::BiCopCDF, VineCopula::BiCopHfunc1),
    lower = -29.99,
    upper = -1.01
  )
)
