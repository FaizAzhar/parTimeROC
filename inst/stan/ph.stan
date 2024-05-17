functions{
  #include /chunks/baselines.stan
  #include /chunks/loglikelihoods.stan
}

data{
  int n;
  int p;
  vector[n] time;
  vector[n] event;
  vector[n] X;
  int baseline;
}

transformed data{
  int is_alpha = 0;
  int is_gamma = 0;
  int is_loc = 0;
  int is_lambda = 0;
  int is_mu = 0;
  int is_sigma = 0;
  int is_omega = 0;
  int is_skew = 0;

  vector[n] zeros;

  if(baseline == 1){ //exponential
    is_lambda = 1;
  }else if(baseline == 2){ // weibull
    is_alpha = 1;
    is_gamma = 1;
  }else if(baseline == 3){ // normal
    is_mu = 1;
    is_sigma = 1;
  }else if(baseline == 4){ // lognormal
    is_mu = 1;
    is_sigma = 1;
  }else if(baseline == 5){ // Gompertz
    is_alpha = 1;
    is_gamma = 1;
  }else if(baseline == 6){ // Skew-Normal
    is_loc = 1;
    is_omega = 1;
    is_skew = 1;
  }
}

  parameters{
    real zbeta;
    array[is_alpha == 0 ? 0 : 1] real<lower=0> alpha;
    array[is_gamma == 0 ? 0 : 1] real<lower=0> gamma;
    array[is_omega ==  0 ? 0 : 1] real<lower=0> omega;
    array[is_lambda == 0 ? 0 : 1] real<lower=0> lambda;
    array[is_mu == 0 ? 0 : 1] real mu;
    array[is_loc == 0 ? 0 : 1] real loc;
    array[is_sigma ==  0 ? 0 : 1] real<lower=0> sigma;
    array[is_skew == 0 ? 0 : 1] real skew;

  }

  model{

    vector[n] lp;
    vector[n] y;
    vector[n] loglik;
    vector[n] lpdf;
    vector[n] lsurv;

    if(p>0){
      zbeta ~ normal(0,10);
      lp = X*zbeta;
    }

    y = time;

    if(baseline == 1){ // exponential
      lambda ~ gamma(1, 1);
      for(i in 1:n){
        lpdf[i] = exponential_lpdf(y[i]|lambda);
        lsurv[i] = exponential_lccdf(y[i]|lambda);
      }
    }else if(baseline == 2){ // Weibull
      alpha ~ gamma(1, 1);
      gamma ~ gamma(1, 1);
      for(i in 1:n){
        lpdf[i] = weibull_lpdf(y[i]|alpha, gamma);
        lsurv[i] = weibull_lccdf(y[i]|alpha, gamma);
      }
    }else if(baseline == 3){ // normal
      mu ~ normal(0, 10);
      sigma ~ gamma(0.1, 1);
      for(i in 1:n){
        lpdf[i] = normal_lpdf(y[i]|mu, sigma);
        lsurv[i] = normal_lccdf(y[i]|mu, sigma);
      }
    }else if(baseline == 4){ // lognormal
      mu ~ gamma(1, 1);
      sigma ~ gamma(0.1, 1);
      for(i in 1:n){
        lpdf[i] = lognormal_lpdf(y[i]|mu, sigma);
        lsurv[i] = lognormal_lccdf(y[i]|mu, sigma);
      }
    }else if(baseline == 5){ // Gompertz
      alpha[1] ~ gamma(1, 1);
      gamma[1] ~ gamma(1, 1);
      for(i in 1:n){
        lpdf[i] = gompertz_lpdf(y[i]|alpha[1], gamma[1]);
        lsurv[i] = gompertz_lccdf(y[i]|alpha[1], gamma[1]);
      }
    }else if(baseline == 6){ // Skew Normal
      loc ~ normal(0,10);
      omega ~ gamma(1, 1);
      skew ~ normal(0,10);
      for(i in 1:n){
        lpdf[i] = skew_normal_lpdf(y[i]|loc, omega, skew);
        lsurv[i] = skew_normal_lccdf(y[i]|loc, omega, skew);
      }
    }


    if(p == 0){
      loglik = event .* lpdf + (1-event) .* lsurv;
    }else{
      loglik = loglik_ph(lpdf, lsurv, event, lp);
    }

    target += sum(loglik);

  }
