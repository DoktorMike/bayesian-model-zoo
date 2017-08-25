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

  vector solowlag(vector x, real lambda, real r, real a) {
    vector[length(x)] y = 0;
    vector[length(x)] b = 0;
    real b0 =lambda*((1-lambda)**r);
    //b=rep(0, length(x))
    b[1]=b0;
    for(i in 1:(length(x)-1)) b[i+1] = b[i]*((r+i-1)/i)*lambda;
    //y=rep(0,length(x));
    for(tt in 1:length(y)) y[tt] = sum(b[1:tt]*x[tt:1]);
    return y;
  }
}
