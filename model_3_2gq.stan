data{
     vector[4993] lnDeltaT_part_1;
     vector[4993] lnTnorm_part_0;
     vector[4993] lnhen_part_0;
     vector[4993] lndrift_part_0;
     vector[4993] lnIM_part_0;
     int N;
     vector[4] alpha;
    array[4993] int ht_part_0;
     vector[4993] lnIM_part_2;
    array[4993] int DS_part_0;
    array[4993] int DS_part_2;
}
parameters{
     vector[max(DS_part_0)] b0;
     real b1;
     real b2;
     ordered[4] cutpoints;
     simplex[4] delta;
}
model{
     vector[4993] eta;
    delta ~ dirichlet( alpha );
    cutpoints ~ normal( 0 , 1.5 );
    b2 ~ normal( 0 , 1 );
    b1 ~ normal( 0 , 1 );
    b0 ~ normal( 0 , 1 );
    for ( i in 1:4993 ) {
        eta[i] = b0[DS_part_0[i]] * lnIM_part_2[i] + b1 * sum(append_row(0, delta)[1:DS_part_0[i]]) + b2 * ht_part_0[i];
    }
    for ( i in 1:4993 ) DS_part_2[i] ~ ordered_logistic( eta[i] , cutpoints );
}
generated quantities{
     vector[N] log_lik;
    for ( i in 1:N ) {
        log_lik[i] = ordered_logistic_lpmf(DS_part_2[i] | b0[DS_part_0[i]] * lnIM_part_2[i] + b1 * sum(append_row(0, delta)[1:DS_part_0[i]]) + b2 * ht_part_0[i], cutpoints);
    }
}


