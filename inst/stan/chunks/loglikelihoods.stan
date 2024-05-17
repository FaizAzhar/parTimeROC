vector loglik_ph(vector lpdf, vector lsurv, vector event, vector lp){
  int n = num_elements(lpdf);
  vector[n] loglik;
  vector[n] lht = lpdf - lsurv;
  loglik = event  .* (lht + lp) +  exp(lp) .* lsurv ;
  return loglik;
}
