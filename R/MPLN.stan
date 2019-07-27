



data{
  int<lower=1> d;      // Dimension of theta
  int<lower=0> N;      //Sample size
  int y[N,d];          //array of Y
  vector[d] mu;
  matrix[d,d] sigma;
  
}
parameters{
  matrix[N,d] theta;

}
model{ //Change values for the priors as appropriate
  for (n in 1:N){
    theta[n,]~multi_normal(mu,sigma);
  }
  for (n in 1:N){
    for (k in 1:d){
      real z;
      z=exp(theta[n,k]);
      y[n,k]~poisson(z);
    }
  }
}


