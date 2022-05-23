functions{
	// flatten_lower_tri: function that returns the lower-tri of a matrix, flattened to a vector
	vector flatten_lower_tri(matrix mat) {
		int n_cols = cols(mat) ;
		int n_uniq = (n_cols * (n_cols - 1)) %/% 2;
		vector[n_uniq] out ;
		int i = 1;
		for(c in 1:(n_cols-1)){
			for(r in (c+1):n_cols){
				out[i] = mat[r,c];
				i += 1;
			}
		}
		return(out) ;
	}
}
data{
	int n ;
	real<lower=0> lkj_prior ;
}
parameters{
	cholesky_factor_corr[n] cfc ;
}
model{
	cfc ~ lkj_corr_cholesky(lkj_prior) ;
}
generated quantities{
	vector[(n*(n-1))%/%2] r = flatten_lower_tri(multiply_lower_tri_self_transpose(cfc)) ;
}
