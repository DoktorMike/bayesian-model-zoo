data {
  int<lower=1> N;
  real y[N];
}

transformed data {
  real D[N];
  real log_unif;
  log_unif = log(N);
  for(i in 1:N) D[i] = (y[i] - mean(y))/sd(y);
}

parameters {
  real mu1;
  real mu2;

  real<lower=0> sigma1;
  real<lower=0> sigma2;
}

// Marginalize out tau and calculate log_p(D | mu1, sd1, mu2, sd2)
transformed parameters {
  vector[N] log_p;
  real mu;
  real sigma;
  log_p = rep_vector(log_unif, N);
  for (tau in 1:N)
    for (i in 1:N) {
      mu = i < tau ? mu1 : mu2;
      sigma = i < tau ? sigma1 : sigma2;
      log_p[tau] = log_p[tau] + normal_lpdf(D[i] | mu, sigma);
    }
}

model {
  mu1 ~ normal(0, 1);
  mu2 ~ normal(0, 1);
  sigma1 ~ normal(0, 1);
  sigma2 ~ normal(0, 1);

  target += log_sum_exp(log_p);
}

generated quantities {
  int<lower=1,upper=N> tau;
  tau = categorical_rng(softmax(log_p));
}
