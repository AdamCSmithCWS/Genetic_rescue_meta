
data {
  int<lower=0> n_obs; // number of observations
  vector[n_obs] He_pre; // vector of pre He values
  vector[n_obs] He_post; // vector of post He values
  vector[n_obs] He_se_pre; // vector of pre SE of He values
  vector[n_obs] He_se_post; // vector of post SE of He values
  
  int<lower=0> n_studies; // number of studies 
  array[n_obs] int<lower=1,upper=n_studies> study; // vector of study indicators



}


parameters {
  real<lower=0> sigma;
  //vector[n_studies] alpha_raw;
  vector[n_studies] beta_raw;
  vector[n_obs] noise_raw;
  //real ALPHA;
  real BETA;
  //real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector[n_obs] true_pre;
  
}

transformed parameters{
  vector[n_studies] beta;
  //vector[n_studies] alpha;
  vector[n_obs] noise;
  vector[n_obs] true_post;
  
  beta = BETA + sigma_beta*beta_raw;
  //alpha = ALPHA + sigma_alpha*alpha_raw;
  noise = sigma*noise_raw;

for(i in 1:n_obs){
  //true_pre[i] = alpha[study[i]];
  true_post[i] = true_pre[i] + beta[study[i]] + noise[i];
}


}

model {
  
  // measurement error
  for(i in 1: n_obs){
   He_pre[i] ~ normal(true_pre[i],He_se_pre[i]);
   He_post[i] ~ normal(true_post[i],He_se_post[i]);
  }
  
  sigma ~ normal(0,1);
  //sigma_alpha ~ normal(0,1);
  sigma_beta ~ normal(0,1);
  
  //alpha_raw ~ normal(0,1);
  beta_raw ~ normal(0,1);
  noise_raw ~ normal(0,1);
  //ALPHA ~ normal(0,1);
  BETA ~ normal(0,1);
  
  }

generated quantities{
  vector[n_obs] log_ratio;
  
  for(i in 1:n_obs){
    vector[2] test;
    test[1] = true_post[i]/true_pre[i];
    test[2] = 500;
    log_ratio[i] = min(test);
    
    test[2] = 1.0/500;
  log_ratio[i] = log(max(test));
}
}
