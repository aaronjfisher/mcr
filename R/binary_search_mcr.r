#Notes - 
# Throughout this code we refer to e_orig as e0, and to e_switch as e1.
# !! = important note to keep track of, or note for future changes.




#' General function for cross-validation
#' 
#' Calculate the cross-validated (CV) loss of a model fitting algorithm \code{fit_fun}.
#'
#' @param y outcome vector
#' @param X covariate matrix
#' @param n_folds number of folds for CV
#' @param fit_fun takes \code{y}, \code{X} as named arguments and returns a prediction model that can be entered into \code{report_fun}.
#' @param report_fun has arguments \code{mod}, in the form of the returned value of \code{fit_fun}, as well as \code{y} and \code{X}, in the same format as the inputs \code{y} and \code{X} to the overall function \code{CV_mod}. \code{report_fun} should return a numeric vector.
#' @param fit_in_sample whether the within-fold performance metrics should also be reported, in addition to the out-of-fold metrics.
#' @param verbose show progress bar
#'
#' @return The mean of \code{report_fun}'s out-of-fold output, averaged across folds. If \code{fit_in_sample=TRUE}, the within-fold output is also shown.
#' @export
CV_mod <- function(y,X, n_folds=min(15,length(y)), fit_fun, report_fun, fit_in_sample=FALSE, verbose=FALSE){

	n_dat <- length(y)
	cv_breaks <- cut(1:n_dat,breaks=n_folds)

	if(verbose) pb <- timerProgressBar(min = 0, max = n_folds)
	for(i in 1:n_folds){
		cv_ind <- cv_breaks == levels(cv_breaks)[i]
		mod_i <- fit_fun(y=y[!cv_ind], X=X[!cv_ind,])
		X_cv_ind <- X[cv_ind,]
		if(class(X_cv_ind) !='matrix') X_cv_ind <- t(X_cv_ind)
		CV_report_i <- report_fun(mod=mod_i, X = X_cv_ind, y = y[cv_ind])
		if(!is.numeric(CV_report_i)) stop('report_fun must return a numeric vector')
		if(i==1){
			in_sample_report <- 
			CV_report <- matrix(NA,n_folds, length(CV_report_i))
			colnames(in_sample_report) <-
			colnames(CV_report) <- names(CV_report_i)
		}
		if(fit_in_sample) in_sample_report[i,] <- report_fun(mod=mod_i, X = X[!cv_ind,], y = y[!cv_ind])
		CV_report[i,] <- CV_report_i

		if(verbose) setTimerProgressBar(pb, i)
	}
	cv_out <- colMeans(CV_report)
	if(!fit_in_sample) return(cv_out)
	if(fit_in_sample) return(c(
		'in_sample'=colMeans(in_sample_report),
		'cv'  = cv_out
	))
	
}










#' Helper function for loop invariant calculations in binary search
#' 
#' For general model classes, returns an augmented ("full") dataset, containing the original dataset, in addition to terms used to approximate or calculate e_switch. Warning: computational resources required are proportional to n times nrep_sample.
#' @param y outcome vector
#' @param X covariate matrix
#' @param n original sample size
#' @param seed seed used for random permuations of the sample
#' @param p1 indeces for the variables to be switched
#' @param nrep_sample setting `nrep_sample=2` corresponds to e_divide to approximate e_switch. Increasing `nrep_sample` further increases the number of terms used in the approximation of e_switch. If `nrep_sample =n,` all permutations are returned.
#' @param warn_dropped whether to give a warning if nrep_sample does not divide evenly into n. In this case, some number of observations (less than nrep_sample) will be dropped.
get_full_sample <- function(y,X,p1,n=length(y),nrep_sample=2,seed=0, warn_dropped=TRUE){
	p <- dim(X)[2]
	p2 <- setdiff(1:p,p1)

	if(nrep_sample<=1) stop('If nrep_sample=1, then no permuted values are created.')
	if(nrep_sample%%1 != 0 | nrep_sample>n) stop('nrep_sample must be an integer <= n.')

	# Steps:
	# 	1 - create several groups of indeces (n_groups of them)
			# note, nrep_sample is also the size of each group
	# 	2 - in each group, compute all possible combinations
	#	3 - tag elements by whether they correspond to permutations or whether they are contained in the original sample `is_perm`.

	combn_inds <- c() # to hold indeces for permuting data.
	n_groups <- floor(n/nrep_sample)
	if(n_groups != n/nrep_sample & warn_dropped) warning(paste(n%%nrep_sample, 'observation(s) dropped before creating permuted data'))
	starts <- 1 + nrep_sample*(0:(n_groups-1)) #beginning of each group
	for(j in 1:n_groups){
		# store all combinations within the group
		group_inds <- which(1:n >= starts[j] & 1:n < starts[j]+nrep_sample)
		combn_inds <- rbind(combn_inds,expand.grid(group_inds,group_inds))
	}
	if(any(duplicated(combn_inds))) stop('dup error')
	# combn_inds
	colnames(combn_inds) <- c('p1','p2')
	

	full_n <- dim(combn_inds)[1]
	if(n>50 & nrep_sample==n){
		if(full_n!=n^2) stop('n error') #workcheck
		warning('full sample is of size n^2 = ',full_n,'. Consider setting "nrep_sample" smaller') 
	}

	old_seed <- rnorm(1, 0, 1000000)
	set.seed(seed)
		pre_mix <- sample(n)
	set.seed(old_seed)
	full_X <- array(NA, dim=c( full_n, dim(X)[2]))
	full_X[,p1] <- (X[pre_mix,])[combn_inds[,1],p1]
	full_X[,p2] <- (X[pre_mix,])[combn_inds[,2],p2]

	return(list(
		is_perm = combn_inds[,1] != combn_inds[,2],
		full_X = full_X,
		full_y = (y[pre_mix])[combn_inds[,2]]
	))
}



#  __  __    _____   _____
# |  \/  |  / ____| |  __ \
# | \  / | | |      | |__) |
# | |\/| | | |      |  _  /
# | |  | | | |____  | | \ \
# |_|  |_|  \_____| |_|  \_\




# Function for estimating MR
#' @param model model for which MR should be evaluated. If applicable, this should match the output of \code{minimize_weighted_loss}, and be an input for \code{get_loss}.
#' @export
#' @rdname mcr_set
get_MR_general <- function(model, precomputed = NA, ...){
	if(is.na(precomputed[1])) precomputed <- precompute_mcr_objects_and_functions(...)
	precomputed$get_MR(model = model, suff_stats=precomputed$suff_stats)
}



#' Compute MR & MCR
#' 
#' Main functions for computing model reliance (MR) and model class reliance (MCR). For further details, see \code{\link{getMCR_internal}}.
#' 
#' @param X a matrix or dataframe of covariates
#' @param y outcome vector
#' @param p1 numeric index marking which columns of X to compute importance for.
#' @param model_class_loss (optional) Signify a preset model class and loss function for which existing code in the package is already written. Currently supports 'linear_mse' for linear regression with the squared error loss; 'kernel_mse' for regression in a reproducing kernel Hilbert space, with the squared error loss; and 'linear_hinge' for linear classification with the hinge loss (y = -1 or 1). If this option is not set, both \code{minimize_weighted_loss} and \code{get_loss} must be set.
#' @param minimize_weighted_loss (optional) a function with named inputs \code{X}, \code{y}, \code{case.weights}, and (optionally) \code{start}. The function should export a model object which minimizes the sum of losses for the dataset (y,X). Any model format is acceptable, as long as it can be read by \code{get_loss}. The (optional) input \code{start} should also take the same format as the output of \code{minimize_weighted_loss}. Regularization constraints should also be accounted for in this function.  
	# !! add examples
#' @param get_loss (optional) a function which takes named inputs \code{model}, \code{y}, and \code{X}, and returns a vector of losses. The input \code{model} should be in the same format of the returned value of \code{minimize_weighted_loss}.
#' @param nrep_sample an integer between 2 and length(y) determining the level of approximation for empirical MR. If nrep_sample=2, e_divide is used. If nrep_sample=length(y), all combinations are computed. See \code{\link{get_full_sample}}.
#' @param loop_ind_args A list of arguments not changing over the binary search (e.g., \code{crossprod(X)}). 
	# (!!) Replace this w/ ... might make things cleaner, but current approach forces all control options to be specified.
#' @return \code{precompute_mcr_objects_and_functions} returns a list containing \itemize{
	#' \item \code{suff_stats} - a list of precomuted objects that do not change over iterations of the MCR binary search loop.
	#' \item \code{get_ld_h_min} - (advanced or internal use) a function that minimizes of equations in the form of  \eqn{h_\gamma}, which has loop-dependent ("ld") arguments. Takes arguments \code{loop_dep_args} (optional), \code{suff_stats} (above), \code{s} (-1 for MCR-, +1 for MCR+), and \code{gam} (coefficient in \eqn{h_\gamma}, see details in main paper). If only computing MR (see get_MR_general), this element may be NULL. 
	#' \item \code{get_MR} - a function taking \code{model} and \code{suff_stats} as input (see above), and which computes the MR of the model using the precomputed object \code{suff_stats}.
#'}
#' @export
#' @rdname mcr_set
precompute_mcr_objects_and_functions <- function(
	# Generic solver input
	X, y, p1,

	# Preset solvers
	model_class_loss = NULL,

	# Generic solvers
	minimize_weighted_loss = NULL,
	get_loss=NULL,
	nrep_sample=2,

	# custom solver
	loop_ind_args=NULL 
	){

	##############################
	##############################
	##############################
	# Determine a function to minimize
	# objectives in the form of "h_{gamma,s}". 
		# (see details in paper)
	

	generic_input_null <- unlist(lapply(
		list(get_loss, minimize_weighted_loss), is.null))

	if( is.null(model_class_loss) &
		is.null(get_loss)){
			stop('either get_loss or model_class_loss must be specified.')
	}

	
	if(is.null(model_class_loss)){
		# minimize h generically
		if(!is.null(loop_ind_args)) warning('overriding loop_ind_args')

		ssts <- get_full_sample(X=X,y=y,p1=p1,nrep_sample=nrep_sample)
		
		get_MR <- function(model, suff_stats){
			e0_value <- mean(get_loss(model = model,
				y=suff_stats$full_y[!suff_stats$is_perm],
				X=suff_stats$full_X[!suff_stats$is_perm,]))
			e1_value <- mean(get_loss(model = model,
				y=suff_stats$full_y[ suff_stats$is_perm],
				X=suff_stats$full_X[ suff_stats$is_perm,]))
			e1_value / e0_value
		}

		get_ld_h_min <- NULL # loop-dependent function for minimizing h
		if(!is.null(minimize_weighted_loss)){
		get_ld_h_min <- function(s, gam, suff_stats, loop_dep_args = NULL){

			#(!!) Could change name of objects "suff_stats" to "io" for loop independent objects.
			ssts <- suff_stats
			weights_gam <- rep(NA,length(ssts$is_perm))
			if(s==-1){ #MCR-
				weights_gam[ ssts$is_perm] <- 1/sum(ssts$is_perm) #avoid using n here, since sample size may != n if e_divide approach is used.
				weights_gam[!ssts$is_perm] <- gam/sum(!ssts$is_perm)
			}else{ #MCR+
				weights_gam[ ssts$is_perm] <- gam/sum(ssts$is_perm)
				weights_gam[!ssts$is_perm] <- 1/sum(!ssts$is_perm)
			}

			if('start' %in% formalArgs(minimize_weighted_loss)){
				h_model <- minimize_weighted_loss(y = ssts$full_y, X=ssts$full_X, case.weights = weights_gam, start = loop_dep_args$mod_lower$h_model) 
					# Ad hoc tests, using mod_lower seems to work better for MCR+
			}else{
				h_model <- minimize_weighted_loss(y = ssts$full_y, X=ssts$full_X, case.weights = weights_gam)
			}
			e0_value <- mean(get_loss(model = h_model,
				y=ssts$full_y[!ssts$is_perm],
				X=ssts$full_X[!ssts$is_perm,]))
			e1_value <- mean(get_loss(model = h_model,
				y=ssts$full_y[ ssts$is_perm],
				X=ssts$full_X[ ssts$is_perm,]))

			coef0 <- 1; coef1 <- 1
			if(s==-1) coef0 <- gam #MCR-
			if(s==1) coef1 <- gam #MCR+

			return(list(
				h_value = coef0 * e0_value + coef1 * e1_value,
				e0_value = e0_value,
				e1_value = e1_value,
				h_model = h_model
				))
		}}


	} else if(model_class_loss == 'linear_mse'){
		RM <- loop_ind_args$reg_matrix

		if( is.null(RM)){
			ssts <- get_suff_stats_lm(y=y,X=X,p1=p1,tol=10^-10)#!! hard coded
			get_ld_h_min <- get_h_min_lm_unregularized
		}
		if(!is.null(RM)){
			ssts <- get_suff_stats_lm(y=y,X=X,p1=p1,tol=10^-10, reg_threshold=Inf, reg_matrix=RM)#!! hard coded
			get_ld_h_min <- get_h_min_lm_regularized
		}

		get_MR <- get_MR_lm

	}else if(model_class_loss == 'kernel_mse'){
		ssts<- get_suff_stats_kernel(y=y, X=X, p1 = p1,
			dat_ref = loop_ind_args$dat_ref,
			kern_fun = loop_ind_args$kern_fun,
			reg_threshold = loop_ind_args$reg_threshold,
			nrep_sample = loop_ind_args$nrep_sample,
			tol = loop_ind_args$tol,
			verbose=loop_ind_args$verbose,
			warn_psd=loop_ind_args$warn_psd,
			warn_duplicate=loop_ind_args$warn_duplicate,
			warn_dropped=loop_ind_args$warn_dropped
			)
		get_ld_h_min <- get_h_min_lm_regularized
		get_MR <- get_MR_lm

	}else if(model_class_loss == 'linear_hinge'){

		return(precompute_mcr_objects_and_functions(
			model_class_loss = NULL,
			X=X, y=y, p1=p1,
			# Generic solver
			minimize_weighted_loss = function(X, y, case.weights, start){

				Xt <- X
				trans_X <- loop_ind_args$transform_x
				if(!is.null(trans_X)) Xt <- t(apply(X, 1, trans_X))
				reg_matrix <- loop_ind_args$reg_matrix
				if(is.null(reg_matrix)) reg_matrix <- diag(dim(Xt)[2])

				optimize_hinge_general(
					X = Xt,
					y = y,
					case.weights = case.weights,
					start = start,
					extra_step_NM = TRUE,
					reg_threshold=loop_ind_args$reg_threshold,
					reg_matrix=reg_matrix)
			},
			get_loss = function(model,y,X){
				Xt <- X
				trans_X <- loop_ind_args$transform_x
				if(!is.null(trans_X)) Xt <- t(apply(X, 1, trans_X))

				hinge_w(w=model, X=Xt, y=y)
			}
		))
	}

	return(list(
		suff_stats = ssts,
		get_ld_h_min = get_ld_h_min,
		get_MR = get_MR
	))
}

#' Calculate MCR- or MCR+
#'
#' Primarily for internal use (see \code{\link{get_empirical_MCR}} for a more user friendly wrapper function).
#' 
#' @param s To search for MCR+, set `s` to 1. To search for MCR-, set `s` to -1. 
#' @param precomputed output from \code{\link{precompute_mcr_objects_and_functions}}
#' @param eps loss threshold on the absolute scale
#' @param tol_mcr tolerance for MCR binary search
#' @param force_lower_0 (TRUE is reccomended) This option can greatly reduce computation time for MCR-, and does not affect MCR+. For MCR-, it tells us to leverage fact that MR values <=1 have a similar interpretation, i.e., null reliance. If MCR- is >=1, this option will have no effect. If MCR- is <=1, this option directs us to not compute to additional precision beyond the fact that it is <=1. 
#' @param warn_lower_0 flags when the search stops early due to finding a well-performing model with MR < 1. As discussed above, this can occur intentionally when force_lower_0=TRUE, and is often not a problem.
#' @param maxiter_BS maximum number of iterations for the binary serach.
#' @param verbose either 'pb' (progress bar),'full', or other for no progress messages.
#' 
#' @return A list containing \itemize{
#' \item \code{mcr} - the MCR value.
#' \item \code{approximation_error_linear} - The difference between the computed bound on MCR, and the most extreme MR value found for a valid model, during the search.
#' \item \code{approximation_error_search} - The largest possible reduction to the approximation error that could be achieved by continuing the search. This is computed based on the fact that our binary search method uses linearizations.
#' \item \code{gam_opt} - (For internal use) final binary search parameter value.
#' \item \code{s} - A value equal to -1 for MCR-, and +1 for MCR. 
#' \item \code{full_model_output} - More detailed description of the final model identified, which will depend on the selected model class.
#' \item \code{searchpath} - Output from steps of the search, used for \code{\link{plot_empirical_MCR_searchpath}}.
#' }
getMCR_internal <- function(s, precomputed, eps,
	tol_mcr = 2^(-20), force_lower_0=TRUE, warn_lower_0=TRUE, maxiter_BS = 400, verbose='pb',...){
	#We search to set gam as low as possible while still meeting cond_move_down


	if(class(eps)!='numeric') stop('eps must be numeric')
	if(class(tol_mcr)!='numeric') stop('tol_mcr must be numeric')
	if(tol_mcr >= 1) stop('tol_mcr must be <1')
	if(!s %in% c(-1,1)) stop('s must be -1 or 1')


	pcof <- precomputed
	if(is.na(pcof[1])) pcof <- precompute_mcr_objects_and_functions(...)

	searchpath<-data.frame('e0'=NA, 'e1'=NA, 'h'=NA, 'gam'=NA)

	get_ld_h_min <- function(gam, add_to_search=TRUE,...){
		mod <- pcof$get_ld_h_min(gam=gam, ...)
		if(add_to_search) searchpath<<-rbind(searchpath,c(
			'e0'=mod$e0_value,
			'e1'=mod$e1_value,
			'h'=mod$h_value,
			'gam'=gam))
		mod
	}
	ssts <- pcof$suff_stats
	if(is.null(get_ld_h_min)) stop('Not enough arguments specified for binary search. Either the set of arguments {get_loss, minimize_weighted_loss}, or model_class_loss must be fully specified.')
	
	# Verbose output
	vcat <- function(...){ if(verbose=='full'){ cat(...) }}

	#Check feasibility
	e0_min_model <- get_ld_h_min(s=1, gam=0, loop_dep_args=NULL,
	 suff_stats=ssts, add_to_search=(s==1))
	if(e0_min_model$e0_value > eps){
		stop('Best model does not meet eps constraint')
	}
	validMR <- e0_min_model$e1_value / e0_min_model$e0_value #baseline to compare against

	#### Conditions
	cond_stop <- function(model_gam){
		((model_gam$h_value == 0 & model_gam$e0_value <= eps) | 
		 (model_gam$h_value >= 0 & model_gam$e0_value == eps))
	}
	cond_move_down <- function(model_gam){
		 (model_gam$h_value  > 0 & model_gam$e0_value <  eps)
	}
	cond_move_up <- function(model_gam){
		 !( cond_move_down(model_gam) | cond_stop(model_gam) )
	}


	#### initialize
	gam_lower <- -1
	gam_upper <- 1
	if(s==1){
		gam_upper <- 0
		mod_upper <- e0_min_model
	}
	if(s==-1){
		gam_lower <- 0 #if force_lower_0, only do convex searches within binary search
		mod_upper <- get_ld_h_min(s=s, gam=gam_upper, suff_stats=ssts)
	}
	mod_lower <- get_ld_h_min(s=s, gam=gam_lower, suff_stats=ssts)


	#### Expand binary search region
	max_expand <- 100
	for(i in 1:max_expand){
		if(!cond_move_up(mod_upper)) break
		gam_upper <- 2*abs(gam_upper+1)
		mod_upper <- get_ld_h_min(s=s, gam=gam_upper, suff_stats=ssts)
	}
	if( (!force_lower_0) | s==1){
		for(i in 1:max_expand){
			if(!cond_move_down(mod_lower)) break
			gam_lower <- -2*abs(gam_lower-1)
			mod_lower <- get_ld_h_min(s=s, gam=gam_lower, suff_stats=ssts)
		}
	}

	gam_opt <- NA
	do_BS_search <- TRUE
	if(cond_move_down(mod_lower)){
		gam_opt <- gam_lower
		model_gam_opt <- mod_lower
		do_BS_search <- FALSE

		if(s==1 | 
			(model_gam_opt$h_value/eps) - gam_opt > 1 | #bound that will be returned
			warn_lower_0) warning('lower limit reached, result may be conservative')
	}
	if(cond_move_up(mod_upper)){
		stop('upper limit reached - size of Rashomon set may be too small? Stopping search to avoid errors.')
	}
	if(gam_lower>gam_upper) stop('search endpoint error')

	#### quick optimality check
	if(cond_stop(mod_lower)){gam_opt <- gam_lower; model_gam_opt <- mod_lower}
	if(cond_stop(mod_upper)){gam_opt <- gam_upper; model_gam_opt <- mod_upper}

	#### Run binary search

	if(verbose == 'pb'){
		cat(c('\n\nMCR-','\n\nMCR+')[(s+3)/2],'binary search progress:\n')
		pb <- txtProgressBar(min=0,max=-log(tol_mcr, base=2))
	}
	
	if(do_BS_search){
		vcat('\n\nStarting binary search\n')
		for(i in 1:maxiter_BS){
			gam_mid <- (gam_lower+gam_upper)/2	
			time_i <- system.time({
				mod_mid <- get_ld_h_min(s=s, gam=gam_mid, 
					loop_dep_args=list(mod_lower =mod_lower, mod_upper=mod_upper), suff_stats=ssts)
			})['elapsed']

			r5 <- function(x) round(x,digits=5)
			vcat('\n sign:',s,'; iter:', i,'; gam:', r5(gam_mid),'; converg:', r5(gam_upper-gam_lower),'; time:',time_i)

			if(cond_stop(mod_mid)){
				gam_upper <- gam_lower <- gam_mid
				mod_upper <- mod_lower <- mod_mid
				break
			}
			if(cond_move_down(mod_mid)){
				gam_upper <- gam_mid
				mod_upper <- mod_mid
				vcat('.. down')
			}
			if(cond_move_up(mod_mid)){
				gam_lower <- gam_mid
				mod_lower <- mod_mid
				vcat('.. up')
			}
			
			if(s==1)  lim_check <- c((mod_upper$h_value/eps - 1) * gam_upper^-1)
			if(s==-1) lim_check <- c((mod_upper$h_value/eps) - gam_upper)
			if(mod_upper$e0_value <= eps) validMR <- mod_upper$e1_value / mod_upper$e0_value


			tol_i <- min(
				abs(gam_upper-gam_lower),
				abs(mod_upper$h_value),
				abs(lim_check-validMR)
				) # Future note - could also use the linear approximation below which captures how much better the binary search we're doing could do, as opposed to how much better *any* method could do. (!!)
			if(verbose=='pb') setTxtProgressBar(pb, -log(max(tol_i,tol_mcr), base=2))
			if(tol_i < abs(tol_mcr)) break
		}
		gam_opt <- gam_upper
		model_gam_opt <- mod_upper
	}


	#### Summarize Results
	vcat('\n\nFound gamma at ',gam_opt,'\n\n')

	if(s==1)  out <- c((model_gam_opt$h_value/eps - 1) * gam_opt^-1)
	if(s==-1) out <- c((model_gam_opt$h_value/eps) - gam_opt)
	if(model_gam_opt$e0_value <= eps){
		validMR <- model_gam_opt$e1_value / model_gam_opt$e0_value
	}else{
		stop('Error in e0 calculation; final value of search should satisfy cond_move_down. Rashomon set may be too small?')
	}
	if(out<0) stop('Negative error')


	#Approximation error compared to a nonlinear approach
	approximation_error_linear <- abs(out - validMR)

	#Approximation error due to binary search resolution, and multiple solutions for get_ld_h_min due to flat parts of the search space.
	approximation_error_search <- NA
	if(mod_lower$h_value > 0 & mod_upper$h_value > 0 & do_BS_search){
		if(sum( c(mod_lower$e0_value, mod_upper$e0_value) > eps) != 1) stop('move condition error') #we should have one model with e0>eps, and one with e0<= eps

		# find the slope and intercept of the line L
		# connecting the coordinates of 
		# mod_upper & mod_lower in e0,e1 space,
		# and take it's intersection with the line e0=eps.
		slope <- (mod_lower$e1_value - mod_upper$e1_value)/(mod_lower$e0_value - mod_upper$e0_value)
		int <- mod_lower$e1_value - slope * mod_lower$e0_value
		worst_case_out <- (int + slope * eps)/eps
		approximation_error_search <- out - worst_case_out
	}


	return(list('mcr'=out,
			'approximation_error_linear'=approximation_error_linear,
			'approximation_error_search'=approximation_error_search,
			'gam_opt'=gam_opt,
			s=s, 
			full_model_output=model_gam_opt,
			searchpath=searchpath[-1,] #drop the NA row
		))

}






# Compute empirical MCR ##documentation set elsewhere.

#' @param eps performance threshold for MCR, on an absolute scale.
#' @param precomputed output from \code{precompute_mcr_objects_and_functions} 
#' @param ... passed to \code{precompute_mcr_objects_and_functions} (and \code{\link{getMCR_internal}} for advanced options).
#' 
#' @return \code{get_empirical_MCR} returns a list containing the MCR range, the epsilon value, and more detailed results (\code{minus} and \code{plus}) for each end of the MCR interval (see \code{\link{getMCR_internal}} for more details).
#' @export
#' @name mcr_set
get_empirical_MCR <- function(eps, precomputed=NA, ...){
	minus <- getMCR_internal(s=-1, eps=eps, precomputed=precomputed,...)
	plus  <- getMCR_internal(s= 1, eps=eps, precomputed=precomputed,...)
	return(list(
		range=c(minus$mcr,plus$mcr),
		minus=minus,
		plus=plus,
		eps=eps
	))
}

#' Plot all models found during MCR binary search, and the associated MR boundaries for any model in the class.
#' 
#' @param emp_mcr output from `get_empirical_MCR`
#' @param eps absolute error value, around which empirical MCR is defined
#' @param add_to_plot whether to add lines to an existing plot, or to create a new plot
#' @param ylim_include vector of numeric values that should be contained in ylim
#' @param xlim_include vector of numeric values that should be contained in xlim
#' @param xlim_1 whether xlim should truncate at value 1 (null reliance)
#' @param show_mcr_range additionally plots brackets corresponding to empirical MCR
#' @param show_all_lines show redundant MR boundaries found in MCR search
#' @param ylab,xlab axis labels, as in \code{\link{plot}}
#' @param fill_col color to shade in region below permorfance threshold. If NA, no shading is applied.
#' @export
plot_empirical_MCR_searchpath <- function(emp_mcr, eps=NA, ylab = NA, xlab = NA, show_all_lines=FALSE, show_mcr_range=FALSE, ylim_include =NULL, add_to_plot =FALSE, xlim_include=NULL, xlim_1 =TRUE, fill_col=NA, ...){

	sp_m <- dplyr::filter(emp_mcr$minus$searchpath, h>=0)
	sp_p <- dplyr::filter(emp_mcr$plus$searchpath, h>=0)
	sp_a <- rbind(sp_m, sp_p)

	###### Setup plot
	sp_a$mr <- sp_a$e1/sp_a$e0
	mr_range <- range(c(1,sp_a$mr, xlim_include))
	xlim <- mr_range
	if(xlim_1) xlim[1] <- max(xlim[1],1)
	ylim_include <- range(c(min(sp_a$e0), eps*1.05, ylim_include), na.rm=TRUE)
	if(is.na(ylab)) ylab <- 'Empirical loss'
	if(is.na(xlab)) xlab <- 'Empirical MR'


	###### Compute boundaries
	mr_global_lim_pts <- c( 
		# points where MR is (conservatively) be globally maximized or minimized.
		# Since the mr_range refers only to valid models with h>=0,
		# these will be (at most) just outside the mr_range (or equal if h=0)
		# we will fill these values later on in with large eps values (see hull_bounds).
		max(-sp_m$gam),
		min(-1/(sp_p$gam[sp_p$gam!=0]))
		)

	mr_seq <- seq(mr_global_lim_pts[1],
				  mr_global_lim_pts[2],length=400)
	mr_seq <- c(mr_seq,sp_a$mr)
	mr_seq <- mr_seq[!duplicated(mr_seq)]
	mr_seq <- mr_seq[order(mr_seq)]

	beyond_lim <-(mr_seq < mr_global_lim_pts[1] | 
				  mr_seq > mr_global_lim_pts[2]  ) 
	if(any(beyond_lim)){
			warning('MR limit error, truncating')
			mr_seq <- mr_seq[!beyond_lim]
	}


	bnds_e0_minus <- apply(sp_m, 1, function(z){
		z['h']/(mr_seq + z['gam'])
	})
	bnds_e0_plus <- apply(sp_p, 1, function(z){
		#includes gam == 0, which is a flat border for e0 min
		z['h']/(z['gam']*mr_seq + 1)
	})


	all_bounds <- cbind(bnds_e0_minus,bnds_e0_plus)
	hull_bounds <- apply(all_bounds,1,max)
	neg_inds <- apply(all_bounds,1,min)<=0
	hull_bounds[neg_inds] <- NA
	hull_bounds[mr_seq %in% mr_global_lim_pts] <-  max(c(1,ylim_include))*50000 #mark points that would be an infinite boundry
	

	###### Plot results
	if(!add_to_plot) plot(c(),
		xlab=xlab, ylab=ylab,
		ylim=ylim_include,
		xlim=xlim, ...)
	if(!is.na(fill_col) & !is.na(eps)){
		poly_ind  <- hull_bounds<=eps
		mr_poly <- mr_seq[poly_ind]
		poly_bounds <- hull_bounds[poly_ind] 
		#ensure that edges are covered if we don't need to search
		#explicitly near eps
		poly_bounds[which(mr_poly==max(mr_poly))] <- eps
		poly_bounds[which(mr_poly==min(mr_poly))] <- eps
		polygon(x=mr_poly, y=poly_bounds,col=fill_col, border=NA)
	}
	points(x=sp_a$mr, y=sp_a$e0, ...)
	abline(v=1,lty=3)
	matlines(x=mr_seq, y=hull_bounds, ...)
	if(show_all_lines){
		all_bounds_NA <- all_bounds
		all_bounds_NA[all_bounds_NA<0] <- NA
		matlines(x=mr_seq, y=all_bounds_NA, lty=3, ...)
	}
	if(show_mcr_range) points(x=emp_mcr$range, y=rep(eps,2), pch=c('[',']'),col='darkgreen',lwd=2)
	if(!is.na(eps)) abline(h=eps,lty=2)


}

# Internal function to plot e_orig versus e_switch, for interactive error checks
# @param emp_mcr output from `get_empirical_MCR`
plot_empirical_MCR_switch_orig <- function(emp_mcr, eps=NA, ylab = NA, xlab = NA, ...){
	sp_m <- dplyr::filter(emp_mcr$minus$searchpath, h>=0)
	sp_p <- dplyr::filter(emp_mcr$plus$searchpath, h>=0)
	sp_a <- rbind(sp_m, sp_p)

	if(is.na(ylab)) ylab <- 'switched empirical loss'
	if(is.na(xlab)) xlab <- 'standard empirical loss'

	plot(x=sp_a$e0,y=sp_a$e1,
		xlab=xlab, ylab=ylab, asp=1, ...)
	if(!is.na(eps)) abline(v=eps)

	#Show boundaries from search lemmas
	abline(b=1,a=0,lty=3)
	abline(b=-emp_mcr$minus$gam_opt,a=emp_mcr$minus$full_model_output$h_value,lty=2)
	abline(b=-1/emp_mcr$plus$gam_opt,a=emp_mcr$plus$full_model_output$h_value/emp_mcr$plus$gam_opt,lty=2)

	legend('topleft',c('MR=1','bounds','eps'),lty=c(3,2,1))
}

# adhoc method to get MCR+ and MCR- from stock optimization methods, and to
# to compare against our binary search result.
#' @import dfoptim
get_bound_tightness <- function(s, eps, get_e0, get_e1, get_norm, reg_threshold, start, threshold_tol = 10^-10, K_start = 50, short_nmk_lim = 20,long_nmk_lim = 5000, maxit_SA = 400){

	invalid_constr <- TRUE
	max_iter <- 20
	iter <- 1
	K <- K_start
	while(invalid_constr & iter < max_iter){
		obj_fun <- function(theta){
			e0_t <- get_e0(theta)
			e1_t <- get_e1(theta)
			norm_t <- get_norm(theta)
			penalty <- 0
			if(e0_t > eps) 				penalty <- penalty + K *(e0_t - eps)
			if(norm_t > reg_threshold) 	penalty <- penalty + K *(norm_t - reg_threshold)

			if(s== -1) return( e1_t/e0_t + penalty)
			if(s==  1) return(-e1_t/e0_t + penalty)
		}

		warm_start <- nmk(par=start, fn=obj_fun,control=list(maxfeval=100))$par

		### we always want SANN to give the same answer, regardless of when it's run in a script
			old_seed <- rnorm(1)*1000000
			set.seed(0) 
				soln1_start_point <- optim(par=warm_start,
					fn=function(proposal){
						nmk(par=proposal, fn=obj_fun,control=list(maxfeval=short_nmk_lim))$value
					},
					method='SANN',
					control=list(maxit=maxit_SA, parscale= abs(warm_start))
				)
			set.seed(old_seed)
		### 

		nmk_soln <- nmk(par=soln1_start_point$par, fn=obj_fun,control=list(maxfeval=long_nmk_lim))
		theta_soln <- nmk_soln$par
		# str(nmk_soln)
		
		iter <- iter+1
		invalid_constr <- (get_norm(theta_soln) > reg_threshold + threshold_tol) |
						  (get_e0(theta_soln)   > eps + threshold_tol)
		K <- max(K,2)*10
	}

	e0_sol <- get_e0(theta_soln)
	e1_sol <- get_e1(theta_soln)
	norm_sol <- get_norm(theta_soln)
	return(list(
			par = theta_soln,
			norm = norm_sol,
			e0 = e0_sol,
			e1 = e1_sol,
			MR = e1_sol / e0_sol,
			full_soln =nmk_soln
		))

}
