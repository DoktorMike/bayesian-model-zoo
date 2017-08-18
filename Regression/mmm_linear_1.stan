
functions {
  real std_normal_lpdf(vector y) {
    return -0.5 * y' * y;
  }

  vector regression_rng(vector beta, matrix x, real sigma) {
    vector[rows(x)] y;
    vector[rows(x)] mu;
    mu = x * beta;
    for (n in 1:rows(x)) y[n] = normal_rng(mu[n], sigma);
    return y;
  }

  vector solowlag(vector x, real lambda, real r) {
    vector[num_elements(x)] y;
    vector[num_elements(x)] b;
    real b0 = lambda*pow(1-lambda, r);
    b[1] = b0;
    for(i in 1:(num_elements(x)-1)) b[i+1] = b[i]*((r+i-1)/i)*lambda;
    for(tt in 1:num_elements(y)) y[tt] = sum(b[1:tt]*x[tt:1]);  // tt:1 range does not WORK!
    //for(tt in 1:num_elements(y)) y[tt] = sum(b[1:tt]*1); // This works!
    return y;
  }
}

data {
  int<lower=0> N;      // Number of data points
  int<lower=0> Kxmi;   // Number of media impressions types
  int<lower=0> Kxmc;   // Number of media clicks types
  int<lower=0> Kxmg;   // Number of media gross types
  int<lower=0> Kxwa;   // Number of weather variables
  vector[N] y;         // The response variable

  matrix[N,7] xswd;    // The weekdays variables
  matrix[N,12] xsm;    // The month variables

  matrix[N,Kxmi] xmi;  // The media variables for impressions
  matrix[N,Kxmc] xmc;  // The media variables for clicks
  matrix[N,Kxmg] xmg;  // The media variables for gross spendings

  matrix[N,Kxwa] xwa;  // All weather variables
}

transformed data {
  vector[N] ynew;
  //ynew = (y - mean(y))/sd(y);
  ynew = y;
}

parameters {
  real<lower=0.01> b0;          // The intercept

  // Seasonality
  vector[7-1] bswd;             // The weekday regression parameters
  vector[12-1] bsm;             // The month regression parameters

  // Media
  vector<lower=0>[Kxmi] bmi;    // The impressions regression parameters
  vector<lower=0>[Kxmc] bmc;    // The clicks regression parameters
  vector<lower=0>[Kxmg] bmg;    // The gross regression parameters

  // Weather
  vector[Kxwa] bwa;             // The weather regression parameters

  real<lower=0> sigma;          // The standard deviation
}

transformed parameters {
  // Declarations
  vector[N] mu;
  vector[7] bswdhat;
  vector[12] bsmhat;

  // The weekday part
  bswdhat[1]=0;
  for (i in 1:(7-1)) bswdhat[i+1] = bswd[i];
  // The month part
  bsmhat[1]=0;
  for (i in 1:(12-1)) bsmhat[i+1] = bsm[i];

  // The mean prediction each timestep
  mu = b0+xswd*bswdhat+xsm*bsmhat+xmi*bmi+xmc*bmc+xmg*bmg+xwa*bwa;
}

model {
  // Seasonality priors
  b0 ~ normal(mean(ynew), sd(ynew));
  bswd ~ normal(0, sd(ynew));
  bsm ~ normal(0, sd(ynew));
  // Media priors
  bmi ~ normal(0, sd(ynew)/sd(xmi));
  bmc ~ normal(0, sd(ynew)/sd(xmc));
  bmg ~ normal(0, sd(ynew)/sd(xmg));
  // Weather priors
  bwa ~ normal(0, sd(ynew)/sd(xwa));
  // Likelihood
  ynew ~ normal(mu, sigma);
}

generated quantities {
  vector[N] yhat;
  vector[N] muhat;
  muhat = b0+xswd*bswdhat+xsm*bsmhat+xmi*bmi+xmc*bmc+xmg*bmg+xwa*bwa;
  for(i in 1:N) yhat[i] = normal_rng(muhat[i], sigma);
  //yhat = yhat*sd(y)+mean(y);
}
