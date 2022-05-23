data{
	int nIds ;
	int nReps ;
	matrix[nIds,nReps] y1 ;
	matrix[nIds,nReps] y2 ;
	matrix[nIds,nReps] y3 ;
}
parameters{
	row_vector[3] means ;
	row_vector<lower=0>[3] sds ;
	vector[nIds] l0_stdnorm ;
	matrix[nIds,3] l123_unique_stdnorm ;
	row_vector<lower=-1,upper=1>[3] r ;
	real<lower=0> noise ;
}
transformed parameters{
	// generate l123_stdnorm such all columns have unit variancce and
	//    cor(l0_stdnorm,l123_stdnorm[,1])==r[1]
	//    cor(l0_stdnorm,l123_stdnorm[,2])==r[2]
	//    cor(l0_stdnorm,l123_stdnorm[,3])==r[3]
	row_vector[3] vars_common = r.^2 ;
	row_vector[3] vars_unique = 1-vars_common ;
	matrix[nIds,3] l123_stdnorm = (
		(
			( rep_matrix(l0_stdnorm,3) .* rep_matrix(r,nIds) )
			+ ( l123_unique_stdnorm .* rep_matrix(sqrt(vars_unique),nIds) )
		)
		./ rep_matrix( sqrt(vars_unique+vars_common) ,nIds)
	) ;
	// scale & shift stdnorms
	matrix[nIds,3] l123 = (
		l123_stdnorm .* rep_matrix(sds,nIds)
		+ rep_matrix(means,nIds)
	);
}
model{
	// prior for latent means/sds & manifest noise
	means ~ std_normal() ; // would want to tweak this to match domain expertise
	sds ~ weibull(2,1) ; // ditto
	noise ~ weibull(2,1) ; // ditto

	// common & unique must be ~ std_normal()
	l0_stdnorm ~ std_normal() ;
	to_vector(l123_unique_stdnorm) ~ std_normal() ;

	// r ~ uniform(-1,1) ; // commented-out bc implied by bounds

	// likelihood
	for(i_nId in 1:nIds){
		y1[i_nId,] ~ normal( l123[i_nId,1] , noise ) ;
		y2[i_nId,] ~ normal( l123[i_nId,2] , noise ) ;
		y3[i_nId,] ~ normal( l123[i_nId,3] , noise ) ;
	}

}

