data{
	int nIds ;
	int nReps ;
	matrix[nIds,nReps] y1 ;
	matrix[nIds,nReps] y2 ;
}
parameters{
	row_vector[2] means ;
	row_vector<lower=0>[2] sds ;
	vector[nIds] common_stdnorm ;
	matrix[nIds,2] unique_stdnorm ;
	real<lower=-1,upper=1> r ;
	real<lower=0> noise ;
}
transformed parameters{
	// generate l1_stdnorm & l2_std_norm such that both have unit variancce and cor(l1,l2)==r
	real w = sqrt(r/(1-r)) ;
	vector[nIds] l1_stdnorm = (common_stdnorm*w + unique_stdnorm[,1]) / sqrt( 1+(w^2) ) ;
	vector[nIds] l2_stdnorm = (common_stdnorm*w + unique_stdnorm[,1]) / sqrt( 1+(w^2) ) ;
	// scale & shift stdnorms
	vector[nIds] l1 = l1_stdnorm*sds[1] + means[1] ;
	vector[nIds] l2 = l2_stdnorm*sds[2] + means[2] ;
}
model{
	// prior for latent means/sds & manifest noise
	means ~ std_normal() ; // would want to tweak this to match domain expertise
	sds ~ weibull(2,1) ; // ditto
	noise ~ weibull(2,1) ; // ditto

	// common & unique must be ~ std_normal()
	common_stdnorm ~ std_normal() ;
	to_vector(unique_stdnorm) ~ std_normal() ;

	// r ~ uniform(-1,1) ; // commented-out bc implied by bounds

	// likelihood
	for(i_nId in 1:nIds){
		y1[i_nId,] ~ normal( l1[i_nId] , noise ) ;
		y2[i_nId,] ~ normal( l2[i_nId] , noise ) ;
	}

}

