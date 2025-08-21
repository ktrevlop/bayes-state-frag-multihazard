data{
  int N;
  array[4] int alpha;
  vector[N] lnTnorm_part_0;
  vector[N] lnhen_part_0;
  vector[N] lndrift_part_0;
  array[N] int DS_part_0;
  vector[N] lnIM_part_0;
  vector[N] lnDeltaT_part_1;
  vector[N] lnIM_part_2;
  array[N] int DS_part_2;
}
parameters{
  real<lower=0> b0;
  real<lower=0> b1;
  ordered[4] cutpoints;
}
model{
  vector[N] eta;
  cutpoints ~ normal( -4 , 2 );
  b1 ~ lognormal( 4 , 2 );
  b0 ~ lognormal( 4 , 2 );
  for ( i in 1:N ) {
    eta[i] = b0 * lnIM_part_2[i] + b1 * exp(lnDeltaT_part_1[i]);
  }
  for ( i in 1:N ) DS_part_2[i] ~ ordered_logistic( eta[i] , cutpoints );
}
generated quantities{
  vector[N] log_lik;
  for ( i in 1:N ) {
    log_lik[i] = ordered_logistic_lpmf(DS_part_2[i] | b0 * lnIM_part_2[i] + b1 * exp(lnDeltaT_part_1[i]), cutpoints);
  }
}
