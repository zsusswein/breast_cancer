//
//
// This stan file fits a logistic growth model

functions {
  real[] logistic_model(real t,       // time
               real[] y,     // system state
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real dXdt[1];
    
    dXdt[1] = (theta[1] * y[1]) * (1 - (y[1] / theta[2]));

    return dXdt;

  }
}
data {
  int<lower = 0> N_obs;           // number of observed measurement times
  int<lower = 0> N_mis;
  real ts[N_obs + N_mis];         // measurement times > 0
  
  int<lower = 0> R;              // number of replicates
  
  real<lower = 0> y0[R];          // initial measured populations
  real<lower = 0> y_obs[N_obs, R];   // measured populations
  
  int t_obs[N_obs];              // observed times
  
  int compute_likelihood;
}
transformed data{
  real x_r[0];
  int x_i[0];
  int<lower = 0> N = N_obs + N_mis;
}
parameters {
  real<lower = 0> theta[2]; // r, K
  real<lower = 0> sigma[1]; // model error
  real<lower = 0> y0_hat[1];    // initial_population
}
transformed parameters{
  real<lower = 0> y_hat[N, 1, R];
  for (r in 1:R){
  y_hat[1:N, :, r] = integrate_ode_rk45(logistic_model, y0_hat, 0, ts, theta, rep_array(0.0, 0), rep_array(0, 0));
  }
}
model {
  theta[{1}] ~ exponential(10);
  theta[{2}] ~ normal(9000, 1000);
  sigma ~ exponential(75);

  y0_hat ~ normal(1000, 500);
  
  if (compute_likelihood == 1){
    
    for (r in 1:R){
    
    y0[r] ~ normal(y0_hat, 500);
    
    for (t in 1:N_obs){
      y_obs[t, r] ~ lognormal(log(y_hat[t_obs[t], 1, r]), sigma); 
    
    }
    }
  }
}
generated quantities{
  real<lower = 0> y_rep[N, 1];
  real<lower = 0> ODE_out[N, 1] = integrate_ode_rk45(logistic_model, y0_hat, 0, ts, theta, rep_array(0.0, 0), rep_array(0, 0));
  
  for (t in 1:N){
    y_rep[t] = lognormal_rng(log(ODE_out[t]), sigma);
  }
}
