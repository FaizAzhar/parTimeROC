

// *****************************************************************
  // Gompertz distribution:

  real gompertz_lpdf(real x, real alpha, real gamma){
    real lpdf = log(alpha) + log(gamma) + gamma*x - alpha*expm1(gamma*x);
    return lpdf;
  }

real gompertz_lccdf(real x, real alpha, real gamma){
  real lsurv = - alpha*expm1(gamma*x);
  return lsurv;
}
