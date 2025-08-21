data{
  int N;
  vector[N] lnDeltaT_part_1;
  vector[N] lnTnorm_part_0;
  vector[N] lnhen_part_0;
  vector[N] lndrift_part_0;
  vector[N] lnIM_part_0;
  vector[4] alpha;
  vector[N] lnIM_part_2;
  array[N] int DS_part_0;
  array[N] int DS_part_2;
}
parameters{
  vector[max(DS_part_0)] b0;
  real b1;
  ordered[4] cutpoints;
  simplex[4] delta;
}
model{
  vector[N] eta;
  delta ~ dirichlet( alpha );
  cutpoints ~ normal( 0 , 1.5 );
  b1 ~ normal( 0 , 1 );
  b0 ~ normal( 0 , 1 );
  for ( i in 1:N ) {
    eta[i] = b0[DS_part_0[i]] * lnIM_part_2[i] + b1 * sum(append_row(0, delta)[1:DS_part_0[i]]);
  }
  for ( i in 1:N ) DS_part_2[i] ~ ordered_logistic( eta[i] , cutpoints );
}
generated quantities{
  vector[N] log_lik;
  for ( i in 1:N ) {
    log_lik[i] = ordered_logistic_lpmf(DS_part_2[i] | b0[DS_part_0[i]] * lnIM_part_2[i] + b1 * sum(append_row(0, delta)[1:DS_part_0[i]]), cutpoints);
  }
}
