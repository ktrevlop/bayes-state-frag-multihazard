data{
  int N;
  vector[N] lnDeltaT_part_1;
  array[N] int DS_part_0;
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
    eta[i] = b0 * lnDeltaT_part_1[i];
  }
  for ( i in 1:N ) DS_part_0[i] ~ ordered_logistic( eta[i] , cutpoints );
}
generated quantities{
  vector[N] log_lik;
  for ( i in 1:N ) {
    log_lik[i] = ordered_logistic_lpmf(DS_part_0[i] | b0 * lnDeltaT_part_1[i], cutpoints);
  }
}
