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
