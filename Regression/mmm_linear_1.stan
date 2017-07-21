data {
  int<lower=0> N;                // Number of data points
  vector[N] y;            // The response variable

  matrix[N,7] xweekday; // The weekdays variables
  matrix[N,12] xmonth;  // The month variables

  matrix[N,3] ximp;     // The media variables
  matrix[N,1] xclicks;  // The media variables
  matrix[N,4] xgross;   // The media variables

  matrix[N,3] xweather; // The weather variables
}

transformed data {
  vector[N] ynew;
  //ynew = (y - mean(y))/sd(y);
  ynew = y;
}

parameters {
  real<lower=0.01> b0;  // The intercept

  // Seasonality
  vector[7-1] bweekday;   // The weekday regression parameters
  vector[12-1] bmonth;    // The month regression parameters

  // Media
  vector<lower=0>[3] bimp;       // The impressions regression parameters
  vector<lower=0>[1] bclicks;    // The clicks regression parameters
  vector<lower=0>[4] bgross;     // The gross regression parameters

  // Weather
  vector[3] bweather;   // The weather regression parameters

  real<lower=0> sigma; // The standard deviation
}

transformed parameters {
  // Declarations
  vector[N] mu;
  vector[7] bweekdayhat;
  vector[12] bmonthhat;

  // The weekday part
  bweekdayhat[1]=0;
  for (i in 1:(7-1)) bweekdayhat[i+1] = bweekday[i];
  // The month part
  bmonthhat[1]=0;
  for (i in 1:(12-1)) bmonthhat[i+1] = bmonth[i];

  // The mean prediction each timestep
  mu = b0+xweekday*bweekdayhat+xmonth*bmonthhat+ximp*bimp+xclicks*bclicks+xgross*bgross+xweather*bweather;
}

model {
  // Priors
  b0 ~ normal(mean(ynew), sd(ynew));
  bweekday ~ normal(0, sd(ynew));
  bmonth ~ normal(0, sd(ynew));
  // Media priors
  bimp ~ normal(0, sd(ynew)/sd(ximp));
  bclicks ~ normal(0, sd(ynew)/sd(xclicks));
  bgross ~ normal(0, sd(ynew)/sd(xgross));
  // Weather priors
  bweather ~ normal(0, sd(ynew)/sd(xweather));
  // Likelihood
  ynew ~ normal(mu, sigma);
}

generated quantities {
  vector[N] yhat;
  yhat = b0+xweekday*bweekdayhat+xmonth*bmonthhat+ximp*bimp+xclicks*bclicks+xgross*bgross+xweather*bweather;
  //yhat = yhat*sd(y)+mean(y);
}
