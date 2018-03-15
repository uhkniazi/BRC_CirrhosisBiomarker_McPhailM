data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    int y[Ntotal]; // response variable binomial distributed
    int ivCoefMap[2]; // mapping variable for the betas for scales
}
// transformed data {
// }
parameters {
  // parameters to estimate in the model
    vector[Ncol] betas; // regression parameters
    //vector<lower=0>[(Ncol-1)] sigmas; // scale parameters
    real<lower=0> tau[2]; // standard deviation for deflections
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  vector[Ncol] betas2;
  betas2[1] = betas[1]; // intercept has independent scale
  betas2[2:ivCoefMap[1]] = betas[2:ivCoefMap[1]] * tau[1];
  betas2[(ivCoefMap[1]+1):ivCoefMap[2]] = betas[(ivCoefMap[1]+1):ivCoefMap[2]] * tau[2];
  mu = X * betas2; 
  mu = inv_logit(mu);
}
model {
  //sigmas ~ cauchy(0, 1);
  tau ~ cauchy(0, 1);
  betas[1] ~ cauchy(0, 10); //prior for the betas
  betas[2:Ncol] ~ normal(0, 1); // 
  
  // likelihood function
  y ~ bernoulli(mu);
}
