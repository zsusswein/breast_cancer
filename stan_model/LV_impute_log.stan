//
// This stan file fits a LV growth model

functions {
  real[] LV_model(real t,       // time
               real[] y,     // system state
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real dXdt[2];
    
    real r_s = x_r[3];
    real K_s = x_r[1];

    real r_r = x_r[4];
    real K_r = x_r[2];

    real alpha = theta[1];
    real beta = theta[2];
    
    dXdt[1] = (r_s * y[1]) * (1 - (y[1] / K_s)- ((alpha * y[2]) / K_s));
    
    dXdt[2] = (r_r * y[2]) * (1 - (y[2] / K_r) - ((beta * y[1]) / K_r));

    return dXdt;

  }
}
data {
  int<lower = 0> N_obs;           // number of observed measurement times
  int<lower = 0> N_mis;           // number of missing measurement times
  real ts[N_obs + N_mis];         // total measurement times > 0

  int<lower = 0> R;              // number of replicates

  real<lower = 0> y0[2, R];   // initial measured populations
  real<lower = 0> y_obs[N_obs, 2, R];    // measured populations
  
  
  int t_obs[N_obs];              // actual observed times
  
  int compute_likelihood;       // boolean switch for prior simulation or true model run
  
  real x_r[4];

}
transformed data{
  int x_i[0];

  int<lower = 0> N = N_obs + N_mis;
}
parameters {
  real theta[2]; // alpha, beta
  real<lower = 0> sigma[1]; // pooled measurement error
  real<lower = 0> y0_hat[2];    // initial_population
}
transformed parameters{
  
  real<lower = 0> y_hat[N, 2, R];
  
  for (r in 1:R){
  y_hat[1:N, 1:2, r] = integrate_ode_rk45(LV_model, y0_hat, 0, ts, theta, x_r, rep_array(0, 0));

  }
}
model {
  
  // measurement error
  sigma ~ exponential(20);
  
  // interaction effect
  theta[{1, 2}] ~ normal(0, 1);
  
  
  for (r in 1:R){
  for (pop in 1:2){
    y0_hat[{pop}] ~ normal(3000, 500);
    
      if (compute_likelihood == 1){
        
        y0[{pop}, r] ~ lognormal(log(y0_hat[{pop}]), sigma);
        for (t in 1:N_obs){
        
          y_obs[t, pop, r] ~ lognormal(log(y_hat[t_obs[t], pop, r]), sigma); 
    
        }
      }
    }
  }
}

generated quantities{
  real y_rep[N, 2, 1];
  real ODE_out[N, 2] = integrate_ode_rk45(LV_model, y0_hat, 0, ts, theta, x_r, rep_array(0, 0));

  for (t in 1:N){
    for (pop in 1:2){
          y_rep[t, pop] = lognormal_rng(log(ODE_out[t, pop]), sigma);
    }
  }
}

