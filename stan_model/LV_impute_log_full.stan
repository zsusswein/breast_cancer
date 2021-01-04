//
// This stan file fits a LV growth model

functions {
  real[] LV_model(real t,       // time
               real[] y,     // system state
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real dXdt[2];
    
    // unpack vars
    real alpha = theta[1];
    real beta = theta[2];
    
    real r_s = theta[3];
    real r_r = theta[4];

    real K_s = theta[5];
    real K_r = theta[6];

    
    dXdt[1] = (r_s * y[1]) * (1 - (y[1] / K_s)- ((alpha * y[2]) / K_s));
    
    dXdt[2] = (r_r * y[2]) * (1 - (y[2] / K_r) - ((beta * y[1]) / K_r));

    return dXdt;

  }
}
data {
  int<lower = 0> N_obs;           // number of observed measurement times
  int<lower = 0> N_mis;           // number of missing measurement times
  real ts[N_obs + N_mis];         // measurement times, indexed from 1

  int<lower = 0> R;              // number of replicates

  real<lower = 0> y0[2, R];   // initial measured populations
  real<lower = 0> y_obs[N_obs, 2, R];    // measured populations
  
  
  int t_obs[N_obs];              // actual observed times
  
  int compute_likelihood;       // boolean switch for prior simulation or true model run
  
  real logistic_pars[8, 1]; // cols: r_s, r_r, K_s, K_r; rows: mu, sigma
  
}
transformed data{
  int x_i[0];
  
  real x_r[0];

  int<lower = 0> N = N_obs + N_mis;
}
parameters {
  //real theta[2]; // alpha, beta
  
  // interaction
  real alpha;
  real beta;
  
  real<lower = 0> sigma[1]; // pooled measurement error
  real<lower = 0> y0_hat[2];    // initial_population
  
  // growth rate
  real<lower = 0> r_s;
  real<lower = 0>r_r;
  
  //carrying capacity
  real<lower = 0> K_s;
  real<lower = 0> K_r;
}
transformed parameters{
  
  real<lower = 0> y_hat[N, 2, R];
  
  // define vector to pass to ODE solver
  real theta[6];
  
  // pack vars into vector
  theta[1] = alpha;
  theta[2] = beta;
  theta[3] = r_s;
  theta[4] = r_r;
  theta[5] = K_s;
  theta[6] = K_r;
  
  for (r in 1:R){
  y_hat[1:N, 1:2, r] = integrate_ode_rk45(LV_model, y0_hat, 0, ts, theta, rep_array(0.0, 0), rep_array(0, 0));

  }
}
model {
  
  // measurement error
  sigma ~ exponential(20);
  
  // interaction effect
  theta[{1, 2}] ~ normal(0, 1);
  
  // growth rate -- from logistic
  r_s ~ normal(logistic_pars[1, 1], logistic_pars[2, 1]);
  r_r ~ normal(logistic_pars[3, 1], logistic_pars[4, 1]);
  
  // carrying capacity -- from logistic
  K_s ~ normal(logistic_pars[5,1], logistic_pars[6, 1]);
  K_r ~ normal(logistic_pars[7, 1], logistic_pars[8, 1]);
  
  // compare observed data to model
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

