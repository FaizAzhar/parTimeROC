#' rate_change
#'
#' @description
#' Function to compute the rate change of the ROC with respect to dx and dt.
#'
#' @param obj A 'fitROC' or 'TimeROC' object.\cr
#' @param t A numeric/vector specifying the time of interest.\cr
#' @param n Number of point on the ROC curve which will be used to check the rate of change.\cr
#' @param type A string specifying the type of analysis used. (Default: 'standard')\cr
#' @param params.x A named vector for the biomarker's parameter.\cr
#' @param params.t A named vector for the time-to-event parameter.\cr
#' @param copula A string indicating the type of copula to be used.\cr
#' @param definition A string indicitaing the definition of ROC to use.\cr
#' @param params.copula A numeric for copula's parameter.\cr
#' @param params.ph A numeric for proportional hazard model's parameter.\cr
#' @param cutoff A numeric to generate total biomarker used.\cr
#' @param seed A numeric to pass in set.seed.\cr
#' @param newx A numeric.\cr
#'
#' @return A list of rate change dt, dx and the angle between these rate of change.
#' @export
#'
#' @examples
#' ## Using 'fitROC' object
#' test <- timeroc_obj("normal-weibull-copula", copula = "gaussian")
#' rr <- rtimeroc(test, n=500,
#'                params.x = c(mean=5, sd=0.8),
#'                params.t = c(shape=1.6, scale=5),
#'               params.copula = -0.3)
#' cc <- timeroc_fit(test, x = rr$x, t = rr$t, event = rr$event)
#' jj <- rate_change(cc, t = c(1,2,10,11))
#'
#' ## Using 'TimeROC' object
#' test <- timeroc_obj("normal-weibull-PH",
#' params.x = c(mean=5, sd=0.8),
#' params.t = c(shape=1.6, scale=5),
#' params.ph = 1)
#' ee <- rate_change(test, t = c(.1,.2))
#'


rate_change <- function(obj, t, n=3, type = 'standard',
                    params.x, params.t, copula, definition = 'c/d',
                    params.copula, params.ph, seed, cutoff = 100, newx){

  ## preprocessing of arguments
  xval <- dat <- cum.haz <- NULL
  x.dist <- t.dist <- iscopula <- NULL

  if(inherits(obj, 'fitTROC')){
    args <- preproc(c(as.list(environment()), call = match.call()),
                    extract_from_fitTROC)
  }else if(inherits(obj, 'TimeROC')){
    args <- preproc(c(as.list(environment()), call = match.call()),
                    TimeROC_predict_process)
  }

  list2env(args, environment())

  tot <- seq(0,1,length.out= n + 2)
  tot <- round(tot * cutoff, 1)
  tot <- tot[-c(1,3+2)]
  x <- xval[tot]

  if(!iscopula){
    dt <- rchange_t(x.dist, cum.haz, params.x, params.t, params.ph, x, t)
    dx <- rchange_x(x.dist, cum.haz, params.x, params.t, params.ph, x, t)


  } else if(iscopula){
    dt <- rchange_t(x.dist, cum.haz, params.x, params.t, params.ph,
                       x, t, iscopula, copula, t.dist, params.copula)
    dx <- rchange_x(x.dist, cum.haz, params.x, params.t, params.ph,
                       x, t, iscopula, copula, t.dist, params.copula)
  }
  angle_dtdx <- calculate_angle(dt,dx)
  res <- list(change_t = dt, change_x = dx, angle = angle_dtdx)

  return(res)
}

# ---------------------------------------------
# Function to calculate the angle of the rate of change ROC
# with respect to x and t
calculate_angle <- function(dt,dx){
  res <- dt[,1]
  dt <- dt[,-1]
  dx <- dx[,-1]
  n <- ncol(dx)
  r <- nrow(dx)
  nama_kol <- c()
  i <- 1
  while(i < n){
    store_room <- c()
    cname <- strsplit(colnames(dt[i]), "_")
    cname[[1]][1] <- 'angle'
    cname <- paste0(cname[[1]][1],"_",cname[[1]][2])
    nama_kol <- append(nama_kol,cname)

    for(j in 1:r){
      b <- c(dt[[i+1]][j],dt[[i]][j]) # ROC point when changing t
      a <- c(dx[[i+1]][j],dx[[i]][j]) # ROC point when changing x

      theta_dx <- atan(a[2]/a[1])
      theta_dt <- atan(b[2]/b[1])

      # radian to degree
      theta_dx <- (theta_dx*180)/pi
      theta_dt <- (theta_dt*180)/pi
      theta <- theta_dx - theta_dt
      # theta <- (theta*180)/pi
      store_room <- rbind(store_room, theta)
    }
    res <- cbind(res,store_room)
    i <- i + 2
  }

  colnames(res) <- c('res',nama_kol)
  return(data.frame(res))
}

# ---------------------------------------------
# Function to calculate the rate change of ROC
# with respect to x
rchange_x <- function(x.dist, cum.haz, params.x, params.t, params.ph,
                         xval, time, iscopula = 0, copula = NULL, t.dist = NULL,
                         params.copula = NULL){

  res <- xval
  if(iscopula == 0){
    for (t in time){
    store_room <- c()
    Ft <- cdf_xt(x.dist$density, cum.haz, params.x, as.vector(params.t),
                 params.ph, xval, t, 2)
    St <- survival_xt(x.dist$density, cum.haz, params.x, as.vector(params.t),
                      params.ph, xval, t, 2)
    for (x in xval){
      F_xt_h <- cdf_xt(x.dist$density, cum.haz, params.x, as.vector(params.t),
                       params.ph, x + 0.0001, t, 1)
      S_xt_h <- survival_xt(x.dist$density, cum.haz, params.x, as.vector(params.t),
                            params.ph, x + 0.0001, t, 1)
      F_xt <- cdf_xt(x.dist$density, cum.haz, params.x, as.vector(params.t),
                     params.ph, x, t, 1)
      S_xt <- survival_xt(x.dist$density, cum.haz, params.x, as.vector(params.t),
                          params.ph, x, t, 1)

      dF <- (F_xt_h - F_xt)/0.0001
      dS <- (S_xt - S_xt_h)/0.0001
      dTP_dx <- -dF/Ft
      dFN_dx <- -dS/St
      store_room <- rbind(store_room,c(dTP_dx, dFN_dx))
      }
    colnames(store_room) <- c(paste0("dTP_",t), paste0("dFN_",t))
    res <- cbind(res,store_room)
    }
  } else if(iscopula == 1){
    for (t in time){
      store_room <- c()
      Ft <- do.call(t.dist$density[[3]], c(list(t), params.t))
      St <- 1 - Ft
      for (x in xval){
        Fx <- do.call(x.dist$density[[3]], c(list(x), params.x))
        fx <- do.call(x.dist$density[[2]], c(list(x), params.x))
        C1 <- copula$density[[4]](Fx, Ft, family = copula$family, par = params.copula)

        dTP_dx <- -(C1 * fx)/Ft
        dFN_dx <- -(fx - C1*fx)/St
        store_room <- rbind(store_room,c(dTP_dx, dFN_dx))
      }
      colnames(store_room) <- c(paste0("dTP_",t), paste0("dFN_",t))
      res <- cbind(res,store_room)
    }
  }
  return(data.frame(res))
}

# ---------------------------------------------
# Function to calculate the rate change of ROC
# with respect to time
rchange_t <- function(x.dist, cum.haz, params.x, params.t, params.ph,
                         xval, time, iscopula = 0, copula = NULL, t.dist = NULL,
                         params.copula = NULL){

  res <- xval
  if(iscopula == 0){
    for (t in time){
      store_room <- c()
      Ft <- cdf_xt(x.dist$density, cum.haz, params.x, as.vector(params.t), params.ph, xval, t, 2)
      St <- survival_xt(x.dist$density, cum.haz, params.x, as.vector(params.t), params.ph, xval, t, 2)
      ft <- integ.Fx(x.dist$density, cum.haz, params.x, as.vector(params.t), params.ph, xval, t, 2)

      for (x in xval){
        F_xt <- cdf_xt(x.dist$density, cum.haz, params.x, as.vector(params.t), params.ph, x, t, 1)
        S_xt <- survival_xt(x.dist$density, cum.haz, params.x, as.vector(params.t), params.ph, x, t, 1)
        integ_1 <- integ.Fx(x.dist$density, cum.haz, params.x, as.vector(params.t), params.ph, x, t, 1)
        integ_2 <- integ.Sx(x.dist$density, cum.haz, params.x, as.vector(params.t), params.ph, x, t, 1)

        dTP_dt <- (F_xt * ft)/(Ft^2) - integ_1/Ft
        dFN_dt <- (S_xt * ft)/(St^2) - integ_2/St
        store_room <- rbind(store_room,c(dTP_dt, dFN_dt))
      }
      colnames(store_room) <- c(paste0("dTP_",t), paste0("dFN_",t))
      res <- cbind(res,store_room)
    }
  } else if (iscopula == 1){
    for (t in time){
      store_room <- c()
      Ft <- do.call(t.dist$density[[3]], c(list(t), params.t))
      St <- 1-Ft
      ft <- do.call(t.dist$density[[2]], c(list(t), params.t))
      for (x in xval){
        Fx <- do.call(x.dist$density[[3]], c(list(x), params.x))
        fx <- do.call(x.dist$density[[2]], c(list(x), params.x))
        C2 <- copula$density[[5]](Fx, Ft, family = copula$family, par = params.copula)
        C <- copula$density[[3]](Fx, Ft, family = copula$family, par = params.copula)
        dTP_dt <- -(Ft * C2 * ft - C * ft)/(Ft^2)
        dFN_dt <- (St * C2 * ft - (Fx - C) * ft)/(St^2)
        store_room <- rbind(store_room,c(dTP_dt, dFN_dt))
      }
      colnames(store_room) <- c(paste0("dTP_",t), paste0("dFN_",t))
      res <- cbind(res,store_room)
    }
  }


  return(data.frame(res))
}

# ---------------------------------------------
# Functions to calculate the F(x,t), S(x,t),
# Pr(T=t,X<x) &  Pr(T=t,X>x) - for PH model
cdf_xt <- function(x.dist, cum.haz, params.x, params.t, params.ph,
                   x, t, fun){
  f <- function(c){
    fx <- do.call(x.dist[[2]], c(list(x = c), params.x))
    Ht <- do.call(cum.haz,c(list(t=t,parms=params.t)))
    return((1-exp(-Ht*exp(c*params.ph)))*fx)
  }
  if(fun == 1){
    return(hcubature(f, lowerLimit = -Inf, upperLimit = x)$integral)
  } else {
    return(hcubature(f, lowerLimit = -Inf, upperLimit = Inf)$integral)
  }

}

integ.Fx <- function(x.dist, cum.haz, params.x, params.t, params.ph,
                     x, t, fun){
  f <- function(c){
    fx <- do.call(x.dist[[2]], c(list(x = c), params.x))
    Ht <- do.call(cum.haz, c(list(t = t,parms=params.t)))
    Ht_h <- do.call(cum.haz, c(list(t = t - 0.0001,parms=params.t)))
    return((exp(-Ht_h*exp(c*params.ph))-exp(-Ht*exp(c*params.ph)))/0.0001 * fx)
  }
  if(fun == 1){
    return(hcubature(f, lowerLimit = -Inf, upperLimit = x)$integral)
  } else {
    return(hcubature(f, lowerLimit = -Inf, upperLimit = Inf)$integral)
  }

}

survival_xt <- function(x.dist, cum.haz, params.x, params.t, params.ph,
                        x, t, fun){
  f <- function(c){
    fx <- do.call(x.dist[[2]], c(list(x = c), params.x))
    Ht <- do.call(cum.haz, c(list(t = t,parms=params.t)))
    return((exp(-Ht*exp(c*params.ph)))*fx)
  }
  if(fun == 1){
    return(hcubature(f, lowerLimit = x, upperLimit = Inf)$integral)
  } else {
    return(hcubature(f, lowerLimit = -Inf, upperLimit = Inf)$integral)
  }

}

integ.Sx <- function(x.dist, cum.haz, params.x, params.t, params.ph,
                     x, t, fun){
  f <- function(c){
    fx <- do.call(x.dist[[2]], c(list(x = c), params.x))
    Ht <- do.call(cum.haz, c(list(t = t,parms=params.t)))
    Ht_h <- do.call(cum.haz, c(list(t = t + 0.0001,parms=params.t)))
    return((exp(-Ht*exp(c*params.ph))-exp(-Ht_h*exp(c*params.ph)))/0.0001 * fx)
  }
  if(fun == 1){
    return(hcubature(f, lowerLimit = x, upperLimit = Inf)$integral)
  } else {
    return(hcubature(f, lowerLimit = -Inf, upperLimit = Inf)$integral)
  }

}
