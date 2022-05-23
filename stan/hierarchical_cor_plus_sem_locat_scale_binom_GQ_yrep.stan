/*
aria: compile = 1
aria: make_local += STAN_NO_RANGE_CHECKS=true
aria: make_local += STAN_CPP_OPTIMS=true
// aria: make_local += CXXFLAGS+=-mtune=native
// aria: make_local += CXXFLAGS+=-march=native
// aria: make_local += STANCFLAGS+=--O1
// aria: make_local += CXXFLAGS+=-O3
// aria: make_local += CXXFLAGS+=-g0
// aria: make_local += STAN_NO_RANGE_CHECKS=true
// aria: make_local += STAN_CPP_OPTIMS=true
// aria: make_local += STANCFLAGS+=--O1
// aria: make_local += STANCFLAGS+=--Oexperimental
// aria: make_local += PRECOMPILED_HEADERS=true
*/

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

	matrix var1_var2_corsem(matrix std_normal_var1, vector cors, matrix std_normal_var2_unique){
		int nrows = rows(std_normal_var2_unique) ;
		int ncols = cols(std_normal_var2_unique) ;
		vector[nrows] vars_common = cors.^2 ;
		vector[nrows] vars_unique = 1-vars_common ;
		matrix[nrows,ncols] std_normal_var2 = (
			(
				(
					std_normal_var1[1:nrows,]
					.* rep_matrix(
						cors
						, ncols
					)
				)
				+ (
					std_normal_var2_unique
					.* rep_matrix(
						sqrt(vars_unique)
						, ncols
					)
				)
			)
			// divide by the square-root of the sum of the squared weights to yield unit-scale variates (since the component variates have unit-scale too)
			./ rep_matrix(
				sqrt( vars_common + vars_unique )
				, ncols
			)
		) ;
		return(std_normal_var2) ;
	}

	matrix shift_and_scale_cols(matrix std_normal_vals, vector shift, vector scale){
		int nrows = rows(std_normal_vals) ;
		int ncols = cols(std_normal_vals) ;
		matrix[nrows,ncols] shifted_and_scaled_vals = (
			rep_matrix(shift,ncols)
			+ (
				std_normal_vals
				.* rep_matrix(scale,ncols)
			)
		) ;
		return(shifted_and_scaled_vals) ;
	}


}

data{
	int<lower=0,upper=1> prior_informed ;
	int<lower=0,upper=1> likelihood_informed ;

	// nI: number of individuals
	int<lower=2> nI ;

	// nXc: number of condition-level predictors
	int<lower=2> nXc ;

	// rXc: number of rows in the condition-level predictor matrix
	int<lower=nXc> rXc ;

	// Xc: condition-level predictor matrix
	matrix[rXc,nXc] Xc ;

	// iXc: which individual is associated with each row in Xc
	array[rXc] int<lower=1,upper=nI> iXc ;

	// nY: num entries in the LA observation vector
	int nY ;

	// T1_Y_gauss: Time1 observations modelled with location-scale Gaussian model
	vector[nY] T1_Y_gauss ;

	// T1_Y_binom: Time1 observations modelled with binomial model
	array[nY] int<lower=0,upper=1> T1_Y_binom ;

	// T2_Y_gauss: Time1 observations modelled with location-scale Gaussian model
	vector[nY] T2_Y_gauss ;

	// T2_Y_binom: Time1 observations modelled with binomial model
	array[nY] int<lower=0,upper=1> T2_Y_binom ;

	// yXc: which row in Xc is associated with each observation in Y
	array[nY] int<lower=1,upper=rXc> yXc ;

}

transformed data{

	matrix[nXc,rXc] Xct = transpose(Xc) ;

}

parameters{

	// real<lower=0> T1_locat_intercept_mean ;
	real T1_locat_intercept_mean ;
	real<lower=0> T1_locat_intercept_sd ;
	vector[nXc-1] T1_locat_coef_mean ;
	vector<lower=0>[nXc-1] T1_locat_coef_sd ;

	matrix[nXc,nI] T1_locat_icoef_indiv_helper ;
	cholesky_factor_corr[nXc] T1_locat_cholfaccorr ;

	// real<lower=0> T1_scale_intercept_mean ;
	real T1_scale_intercept_mean ;
	real<lower=0> T1_scale_intercept_sd ;
	vector[nXc-1] T1_scale_coef_mean ;
	vector<lower=0>[nXc-1] T1_scale_coef_sd ;

	vector<lower=-1,upper=1>[nXc] T1_locat_scale_cors ;
	matrix[nXc,nI] T1_scale_icoef_indiv_unique ;

	// real<lower=0> T1_binom_intercept_mean ;
	real T1_binom_intercept_mean ;
	real<lower=0> T1_binom_intercept_sd ;
	vector[nXc-1] T1_binom_coef_mean ;
	vector<lower=0>[nXc-1] T1_binom_coef_sd ;

	vector<lower=-1,upper=1>[nXc] T1_locat_binom_cors ;
	matrix[nXc,nI] T1_binom_icoef_indiv_unique ;


	// real<lower=0> T2_locat_intercept_mean ;
	real T2_locat_intercept_mean ;
	real<lower=0> T2_locat_intercept_sd ;
	vector[nXc-1] T2_locat_coef_mean ;
	vector<lower=0>[nXc-1] T2_locat_coef_sd ;

	vector<lower=-1,upper=1>[nXc] T1_T2_locat_cors ;
	matrix[nXc,nI] T2_locat_icoef_indiv_unique ;

	// real<lower=0> T2_scale_intercept_mean ;
	real T2_scale_intercept_mean ;
	real<lower=0> T2_scale_intercept_sd ;
	vector[nXc-1] T2_scale_coef_mean ;
	vector<lower=0>[nXc-1] T2_scale_coef_sd ;

	vector<lower=-1,upper=1>[nXc] T1_T2_scale_cors ;
	matrix[nXc,nI] T2_scale_icoef_indiv_unique ;

	// real<lower=0> T2_binom_intercept_mean ;
	real T2_binom_intercept_mean ;
	real<lower=0> T2_binom_intercept_sd ;
	vector[nXc-1] T2_binom_coef_mean ;
	vector<lower=0>[nXc-1] T2_binom_coef_sd ;

	vector<lower=-1,upper=1>[nXc] T1_T2_binom_cors ;
	matrix[nXc,nI] T2_binom_icoef_indiv_unique ;

}
generated quantities{
	// Y_gauss_rep: posterior-predictive LA observations
	vector[nY] T1_Y_gauss_rep ;
	array[nY] int T1_Y_binom_rep ;
	vector[nY] T2_Y_gauss_rep ;
	array[nY] int T2_Y_binom_rep ;

	{ // local environment to avoid saving intermediate quantities
		// corStdNorms from cors & std-normal helpers ----
		matrix[nXc,nI] T1_locat_icoef_indiv_corStdNorms = (
			T1_locat_cholfaccorr
			* T1_locat_icoef_indiv_helper
		) ;

		// SEMs ----
		matrix[nXc,nI] T1_scale_icoef_indiv_corStdNorms = var1_var2_corsem(
			T1_locat_icoef_indiv_corStdNorms // std_normal_var1
			, T1_locat_scale_cors // cors
			, T1_scale_icoef_indiv_unique // std_normal_var2_unique
		) ;
		matrix[nXc,nI] T1_binom_icoef_indiv_corStdNorms = var1_var2_corsem(
			T1_locat_icoef_indiv_corStdNorms // std_normal_var1
			, T1_locat_binom_cors // cors
			, T1_binom_icoef_indiv_unique // std_normal_var2_unique
		) ;
		matrix[nXc,nI] T2_locat_icoef_indiv_corStdNorms = var1_var2_corsem(
			T1_locat_icoef_indiv_corStdNorms // std_normal_var1
			, T1_T2_locat_cors // cors
			, T2_locat_icoef_indiv_unique // std_normal_var2_unique
		) ;
		matrix[nXc,nI] T2_scale_icoef_indiv_corStdNorms = var1_var2_corsem(
			T1_scale_icoef_indiv_corStdNorms // std_normal_var1
			, T1_T2_scale_cors // cors
			, T2_scale_icoef_indiv_unique // std_normal_var2_unique
		) ;
		matrix[nXc,nI] T2_binom_icoef_indiv_corStdNorms = var1_var2_corsem(
			T1_binom_icoef_indiv_corStdNorms // std_normal_var1
			, T1_T2_binom_cors // cors
			, T2_binom_icoef_indiv_unique // std_normal_var2_unique
		) ;


		// Shifting & scaling ----
		matrix[nXc,nI] T1_locat_icoef_indiv = shift_and_scale_cols(
			T1_locat_icoef_indiv_corStdNorms // std_normal_vals
			, append_row( T1_locat_intercept_mean , T1_locat_coef_mean ) // shift
			, append_row( T1_locat_intercept_sd , T1_locat_coef_sd ) // scale
		) ;
		matrix[nXc,nI] T1_scale_icoef_indiv = shift_and_scale_cols(
			T1_scale_icoef_indiv_corStdNorms // std_normal_vals
			, append_row( T1_scale_intercept_mean , T1_scale_coef_mean ) // shift
			, append_row( T1_scale_intercept_sd , T1_scale_coef_sd ) // scale
		) ;
		matrix[nXc,nI] T1_binom_icoef_indiv = shift_and_scale_cols(
			T1_binom_icoef_indiv_corStdNorms // std_normal_vals
			, append_row( T1_binom_intercept_mean , T1_binom_coef_mean ) // shift
			, append_row( T1_binom_intercept_sd , T1_binom_coef_sd ) // binom
		) ;
		matrix[nXc,nI] T2_locat_icoef_indiv = shift_and_scale_cols(
			T2_locat_icoef_indiv_corStdNorms // std_normal_vals
			, append_row( T2_locat_intercept_mean , T2_locat_coef_mean ) // shift
			, append_row( T2_locat_intercept_sd , T2_locat_coef_sd ) // scale
		) ;
		matrix[nXc,nI] T2_scale_icoef_indiv = shift_and_scale_cols(
			T2_scale_icoef_indiv_corStdNorms // std_normal_vals
			, append_row( T2_scale_intercept_mean , T2_scale_coef_mean ) // shift
			, append_row( T2_scale_intercept_sd , T2_scale_coef_sd ) // scale
		) ;
		matrix[nXc,nI] T2_binom_icoef_indiv = shift_and_scale_cols(
			T2_binom_icoef_indiv_corStdNorms // std_normal_vals
			, append_row( T2_binom_intercept_mean , T2_binom_coef_mean ) // shift
			, append_row( T2_binom_intercept_sd , T2_binom_coef_sd ) // binom
		) ;

		// dot products ----
		row_vector[rXc] T1_locat_cond = columns_dot_product(
			T1_locat_icoef_indiv[,iXc]
			, Xct
		) ;
		row_vector[rXc] T1_scale_cond = sqrt(exp(columns_dot_product(
			T1_scale_icoef_indiv[,iXc]
			, Xct
		))) ;
		row_vector[rXc] T1_binom_cond = columns_dot_product(
			T1_binom_icoef_indiv[,iXc]
			, Xct
		) ;
		row_vector[rXc] T2_locat_cond = columns_dot_product(
			T2_locat_icoef_indiv[,iXc]
			, Xct
		) ;
		row_vector[rXc] T2_scale_cond = sqrt(exp(columns_dot_product(
			T2_scale_icoef_indiv[,iXc]
			, Xct
		))) ;
		row_vector[rXc] T2_binom_cond = columns_dot_product(
			T2_binom_icoef_indiv[,iXc]
			, Xct
		) ;

		// Likelihood ----
		for(i_nY in 1:nY){
			T1_Y_gauss_rep[i_nY] = normal_rng(
				T1_locat_cond[yXc[i_nY]]
				, T1_scale_cond[yXc[i_nY]]
			) ;
			T1_Y_binom_rep[i_nY] = bernoulli_logit_rng(
				T1_binom_cond[yXc[i_nY]]
			) ;
			T2_Y_gauss_rep[i_nY] = normal_rng(
				T2_locat_cond[yXc[i_nY]]
				, T2_scale_cond[yXc[i_nY]]
			) ;
			T2_Y_binom_rep[i_nY] = bernoulli_logit_rng(
				T2_binom_cond[yXc[i_nY]]
			) ;
		}

	}
}
