set_baseline <- function(baseline){
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
