# preamble ----

# install dependencies & source helpers
source('_imports.r')

# load tidyverse and a magrittr pipe
library(tidyverse)
`%<>%` = magrittr::`%<>%`

#make RStudio use cmdstan for syntax checking
enable_rstudio_stan_syntax_check()

#start parallel compile jobs in the background
compile_stan_files_as_parallel_jobs(path='stan')

# generate the data structure ----
nI = 30
reps = 30

(
	expand_grid(
		id = 1:nI
		, ab = factor(c('A','B'))
		, rep = 1:reps
	)
	%>% mutate(
		t1_y_gauss = 0 #just a dummy value for now; we'll generate data using Stan
		, t1_y_binom = 0 #ditto
		, t2_y_gauss = 0 #ditto
		, t2_y_binom = 0 #ditto
	)
) ->
	dat

# add contrasts
(
	dat
	%>% select(-contains('_y_'),-rep)
	%>% distinct()
	%>% mutate(
		contrasts = get_contrast_matrix_rows_as_list(
			data = .
			, formula = ~ ab
			, contrast_kind = halfsum_contrasts
		)
	)
) ->
	Xc_with_vars

#quick view of the contrasts
(
	Xc_with_vars
	%>% unnest(contrasts)
	%>% select(-id)
	%>% distinct()
	%>% View()
)

#package the data for Stan
data_for_stan = lst(

	# Xc: condition-level predictor matrix
	# matrix[rXc,nXc] Xc ;
	Xc = (
		Xc_with_vars
		%>% select(contrasts)
		%>% unnest(contrasts)
		%>% as.matrix()
	)

	# iXc: which individual is associated with each row in Xc
	# array[rXc] int<lower=1,upper=nI> iXc ;
	, iXc = (
		Xc_with_vars
		%>% pull(id)
		%>% as.factor()
		%>% as.numeric()
	)

	# T1_Y_gauss: observations modelled with location-scale Gaussian model
	# vector[nY] T1_Y_gauss ;
	, T1_Y_gauss = dat$t1_y_gauss

	# T1_Y_binom: observations modelled with binomial model
	# array[nY] int<lower=0,upper=1> T1_Y_binom ;
	, T1_Y_binom = dat$t1_y_binom

	# T2_Y_gauss: observations modelled with location-scale Gaussian model
	# vector[nY] T2_Y_gauss ;
	, T2_Y_gauss = dat$t2_y_gauss

	# T2_Y_binom: observations modelled with binomial model
	# array[nY] int<lower=0,upper=1> T2_Y_binom ;
	, T2_Y_binom = dat$t2_y_binom

	# yXc: which row in Xc is associated with each observation in Y
	# array[nY] int<lower=1,upper=rXc> yXc ;
	, yXc = (
		Xc_with_vars
		%>% mutate(Xc_row = 1:n())
		# right-join with dat to preserve dat's row order
		%>% right_join((
			dat
			%>% mutate(
				dat_row = 1:n()
			)
		))
		%>% arrange(dat_row)
		%>% pull(Xc_row)
	)


	# nI: number of individuals
	# int<lower=2> nI ;
	, nI = length(unique(dat$id))

	# nXc: number of condition-level predictors
	# int<lower=2> nXc ;
	, nXc = ncol(Xc)

	# rXc: number of rows in the condition-level predictor matrix Xc
	# int<lower=nXc> rXc ;
	, rXc = nrow(Xc)

	# nY: num entries in the LA observation vector
	# int nY ;
	, nY = length(T1_Y_gauss)

)


# Sample the prior ----
data_for_stan$prior_informed = 1
data_for_stan$likelihood_informed = 0

# load the sampling model (make sure compile jobs were successful!)
mod = cmdstanr::cmdstan_model('stan/hierarchical_cor_plus_sem_locat_scale_binom.stan')

prior_predictive_output = sample_mod(
	data = data_for_stan
	, mod = mod
	, max_treedepth = 11 # 10 is default
	, refresh_perc = 10
	, init = 2 # default of 2 is good for Gaussian likelihoods, lower may be necessary for binomial
)

#check diagnostics (added as an attr; really should subclass instead)
prior_predictive_output %<>% add_draws_and_diagnostics_attr()
print(attr(prior_predictive_output,'dd')$sampler_diagnostics_across_chain_summary)
(
	attr(prior_predictive_output,'dd')$par_summary
	%>% select(rhat,contains('ess'))
	%>% summary()
)

# Generate yreps from a single prior draw ----
gq_mod = cmdstanr::cmdstan_model('stan/hierarchical_cor_plus_sem_locat_scale_binom_GQ_yrep.stan')

(
	attr(prior_predictive_output,'dd')$draws
	%>% filter(.draw==1)
	%>% posterior::as_draws_array()
) ->
	prior_predictive_draw_for_yrep

prior_predictive_draw_gq = gq_mod$generate_quantities(
	data = data_for_stan
	, fitted_params = prior_predictive_draw_for_yrep
)

(
	prior_predictive_draw_gq$draws(format='draws_df')
	%>% as_tibble()
	%>% select(
		.chain,.iteration,.draw
		, contains('Y_gauss_rep')
		, contains('Y_binom_rep')
	)
	%>% pivot_longer(
		cols = -c(.chain,.iteration,.draw)
		, names_to = 'variable'
	)
	%>% separate(
		variable
		, into=c('variable','index','dummy')
		, sep=c('[\\[\\]]')
		, fill = 'right'
		, convert = TRUE
	)
	%>% select(-dummy,-.chain,-.iteration,-.draw)
	%>% arrange(variable,index)
	%>% pivot_wider(names_from=variable)
) ->
	prior_predictive_draw_yreps

data_for_stan$T1_Y_gauss = prior_predictive_draw_yreps$T1_Y_gauss_rep
data_for_stan$T1_Y_binom = prior_predictive_draw_yreps$T1_Y_binom_rep
data_for_stan$T2_Y_gauss = prior_predictive_draw_yreps$T2_Y_gauss_rep
data_for_stan$T2_Y_binom = prior_predictive_draw_yreps$T2_Y_binom_rep

# sample given the yrep of the prior draw ----
data_for_stan$likelihood_informed = 1

post = sample_mod(
	data = data_for_stan
	, mod = mod
	, max_treedepth = 11 # 10 is default
	, refresh_perc = 10
	, init = 2
)

#check diagnostics
post %<>% add_draws_and_diagnostics_attr()
print(attr(post,'dd')$sampler_diagnostics_across_chain_summary)
(
	attr(post,'dd')$par_summary
	%>% select(rhat,contains('ess'))
	%>% summary()
)

# plot_par(
# 	par=c('locat_intercept_mean','locat_coef_mean')
# 	, true = prior_predictive_draw_for_yrep
# )

plot_par(
	post
	, par_subst = '_mean'
	, true = prior_predictive_draw_for_yrep
)

plot_par(
	post
	, par_subst = '_sd'
	, true = prior_predictive_draw_for_yrep
)


# plot_par(
# 	par = 'locat_cors'
# 	, true = prior_predictive_draw_for_yrep
# )
plot_par(
	post
	, par_substr = 'cors'
	, true = prior_predictive_draw_for_yrep
)

