set_baseline <- function(baseline){
  if(!is.character(baseline)) stop("Incorrect baseline")
  baseline <- switch(baseline,
                     "exponential" = 1,
                     "weibull" = 2,
                     "gaussian" = 3,
                     "normal" = 3,
                     "lognormal" = 4,
                     "gompertz" = 5,
                     "skewnormal" = 6
  )
  return(baseline)
}

setup.pch <- function(t, breakpoints){
  # breakpoints <- c(0,breakpoints,Inf)
  pos <- findInterval(t, vec = breakpoints)
  exposure_time <- numeric(0)
  # Create event matrix
  mat_event <- matrix(0, nrow = length(t), ncol = (length(breakpoints)-1))
  for(i in seq_along(t)){
    mat_event[i,pos[i]] <- 1
    exposure_time[i] <- t[i] - breakpoints[pos[i]]
  }

  # Create time exposure matrix
  mat_exposure <- matrix(1:(length(breakpoints)-1), nrow = length(t),
                         ncol = (length(breakpoints)-1), byrow = TRUE)
  full_interval <- diff(breakpoints)
  for(i in seq_along(t)){
    mat_exposure[i,] <- mat_exposure[i,] - pos[i]
  }

  mat_exposure[mat_exposure > 0] <- 0
  mat_exposure[mat_exposure < 0] <- full_interval[col(mat_exposure)][mat_exposure < 0]
  for(i in seq_along(t)){
    mat_exposure[i,pos[i]] <- exposure_time[i]
  }

  return(list(mat_event = mat_event, mat_exposure = mat_exposure))
}
