data{
  int N;
  vector[N] lnIM;
  array[N] int DS;
}
parameters{
  real b0;
  ordered[4] cutpoints;
}
model{
  vector[N] eta;
  cutpoints ~ normal( 0 , 1.5 );
  b0 ~ normal( 0 , 1 );
  for ( i in 1:N ) {
    eta[i] = b0 * lnIM[i];
  }
  for ( i in 1:N ) DS[i] ~ ordered_probit( eta[i] , cutpoints );
}
generated quantities{
  vector[N] log_lik;
  for ( i in 1:N ) {
    log_lik[i] = ordered_probit_lpmf(DS[i] | b0 * lnIM[i], cutpoints);
  }
}
