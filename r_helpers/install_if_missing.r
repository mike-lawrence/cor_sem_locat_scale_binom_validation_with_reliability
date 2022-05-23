cat('install_if_missing() now available as a function\n')
#' Installs any packages not already installed
#' @examples
#' \dontrun{
#' install_if_missing(c('tidyverse','github.com/stan-dev/cmdstanr'))
#' }
install_if_missing = function(pkgs){
	missing_pkgs = NULL

	for(this_pkg in pkgs){
		path = NULL
		try(
			path <-	find.package(basename(this_pkg),quiet=T,verbose=F)
			, silent = T
		)
		if(length(path)==0){
			missing_pkgs = c(missing_pkgs,this_pkg)
		}
	}
	cran_missing = missing_pkgs[!grepl('github.com/',fixed=T,missing_pkgs)]
	if(length(cran_missing)>0){
		message('The following required but uninstalled CRAN packages will now be installed:\n',paste(cran_missing,collapse='\n'))
		install.packages(cran_missing)
	}
	github_missing = missing_pkgs[grepl('github.com/',fixed=T,missing_pkgs)]
	github_missing = gsub('github.com/','',github_missing)
	if(length(github_missing)>0){
		message('The following required but uninstalled Github packages will now be installed:\n',paste(github_missing,collapse='\n'))
		remotes::install_github(github_missing)
	}
	invisible()
}
