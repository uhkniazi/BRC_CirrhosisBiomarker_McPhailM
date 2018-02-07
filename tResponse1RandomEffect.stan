data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each observation to group 1 
  real y[Ntotal]; // response variable normally distributed
  // additional parameters
  real gammaShape; // hyperparameters for the gamma distribution 
  real gammaRate;
  real intercept;
  real intercept_sd;
  real<lower=1> nu;
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  real betas; // regression intercept
  real<lower=0> sigmaRan1; // random effect standard deviation for group 1
  real<lower=0> sigmaPop; // population standard deviation
  vector[Nclusters1] rGroupsJitter1; // number of random jitters for each level of cluster/group 1
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = betas + rGroupsJitter1[NgroupMap1];
}
model {
  //nu ~ exponential(1/29.0);
  sigmaRan1 ~ gamma(gammaShape, gammaRate);
  sigmaPop ~ gamma(gammaShape, gammaRate);
  betas ~ normal(intercept, intercept_sd);
  // random effects sample
  rGroupsJitter1 ~ normal(0, sigmaRan1);
  // likelihood function
  y ~ student_t(nu, mu, sigmaPop);
}
