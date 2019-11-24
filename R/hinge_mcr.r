#Notes - 
# Throughout this code we refer to e_orig as e0, and to e_switch as e1.
# !! = important note to keep track of, or note for future changes.



###### Weighted loss functions, objective, gradient

#' @name hinge_set
hinge_pred <- function(pred, y, m=1){
	(m-pred*y)*(m-pred*y>0)
}

#' Calculate the hinge loss
#' 
#' \code{hinge_pred} computes the in-sample hinge loss, given a set of predictions. \code{hinge_w} computes predictions, and passes these to \code{hinge_pred} to compute the in-sample hinge loss.
#' 
#' @param w weight vector (linear classification model). First element should be an interecept term. For example, the output of \code{\link{optimize_hinge_general}}.
#' @param X covariate matrix, should not include a constant (intercept) column.
#' @param y outcome vector, with elements -1 or 1
#' @param m margin parameter
#' @param ... passed to hinge_pred
#' @param pred y predictions from the model w
#' @export
#' @name hinge_set
hinge_w <- function(w, X, ...){
	if(dim(X)[2]+1 != length(w)) stop(paste('model is of dimension',length(w), 'but data is of dimension', dim(X)[1],'by',dim(X)[2]))

	hinge_pred(pred=w[1]+ X %*% w[-1],...)
}
get_hinge_grad <- function(X,y,w, reg_threshold, reg_matrix, K, case.weights){
	pos_ind <- c(hinge_w(w=w, X=X, y=y, m=1) > 0)
	#1/n sum (m-yi xi w), deriv is
	#1/n sum -yi xi where xi is a vector
	penalty_grad <- rep(0,length(w))
	norm_w <- t(w[-1])%*%reg_matrix %*% w[-1]
	if(norm_w > reg_threshold) penalty_grad[-1] <- K * (t(reg_matrix) + reg_matrix) %*% w[-1]

	n_y_pi_cw <- -y * pos_ind * case.weights 
	err_pos <- c(sum(n_y_pi_cw), crossprod(n_y_pi_cw, X)) #first entry is for the intercept

	err_pos + penalty_grad
}
get_hinge_obj <- function(X,y,w, reg_threshold, reg_matrix, K, case.weights){
	un_penalized <- sum(
			case.weights * 
			c(hinge_w(w=w, X=X, y=y, m=1))
		)
	norm_w <- w[-1] %*% reg_matrix %*% w[-1]
	penalty <- 0
	if(norm_w>reg_threshold) penalty <- K* (norm_w-reg_threshold)
	un_penalized + penalty
}


######### Workcheck
	# library(numDeriv)

	# p <- 15
	# N <- 50
	# X <- matrix(rnorm(N*p),N,p)
	# B_gen <- c(1,(1:p)/p)
	# logity <- c(cbind(1,X) %*% B_gen)
	# y <- rbinom(N,1,prob = 1/(1+exp(-logity)))
	# K <- 100
	# w_eval <- c(0,(1:p)/p)
	# reg_matrix <- matrix(log(1+1:p^2), p, p )
	# crossprod(w_eval)

	# for(test_cw in 1:2){
	# for(test_reg in 1:2){
	# for(test_w in 1:2){
		
	# 	if(test_cw==1) cw <- rnorm(N)
	# 	if(test_cw==2) cw <- runif(N)
		
	# 	if(test_w==1) w_eval <- c(0,(1:p)/p)
	# 	if(test_w==2) w_eval <- B_gen
		
	# 	if(test_reg==1) regt <- Inf
	# 	if(test_reg==2) regt <- crossprod(w_eval)*.75

	# 	fn_K <- function(w) get_hinge_obj(X=X,y=y,w=w, reg_threshold=regt, reg_matrix=reg_matrix, K=K, case.weights=cw)
	# 	gr_K <- function(w) get_hinge_grad(X=X,y=y,w=w, reg_threshold=regt, reg_matrix=reg_matrix, K=K, case.weights=cw)
	# 	fn_K(w_eval)
	# 	# gr_K(w_eval)

	# 	diffs <- abs(grad(fn_K,w_eval) - gr_K(w_eval))
	# 	str(max(abs(diffs)) / max(abs(gr_K(w_eval))))
	# }}}
#########




#' Optimize a linear classifier with possibly negative observation weights
#'
#'  
#' The optimization uses simulated annealing (SA), with a gradient-based search at each step. (Used in our paper, but not the primary focus of this package; the default option of ignore_below_zero=TRUE may cause problems in uses of this function beyond the scope of MCR, outside of this package.)
#' 
#' @param y outcome vector with elements -1 or 1.
#' @param X covariate matrix (n x p), which should not contain an interecept or constant column.
#' @param reg_matrix Matrix R, where w'Rw is the penalty of a coefficient vector w.
#' @param n sample size
#' @param case.weights vector of numeric multipliers (possibly negative) for each observation, when computing loss summation.
#' @param reg_threshold Value for w'Rw above which the penalty is applied
#' @param K penalty multiplier for w'Rw
#' @param p covariate dimension
#' @param start starting value for coefficient vector
#' @param constr_tol buffer for regularization, after which a penalty is applied
#' @param sparse_warning warn user if intercept column is detected
#' @param maxit_SA number of iterations for SA steps
#' @param extra_step_NM whether to follow the SA search with an additional Nelder-Mead search
#' @param short_cg_lim tuning parameter used to set initial value of the search
#' @param long_cg_lim how many gradient-based steps to take at each SA iteration
#' @param ignore_below_zero stop SA search if a value is discovered below zero. This is irrelevant if weights are positive, and is useful within the MCR binary search, but may cause problems for other applications of this function beyond the computation of MCR.
#' 
#' @import optimr dfoptim stats
#' @return a linear coefficient vector
#'
#' @export
optimize_hinge_general <- function(y, X, reg_matrix=diag(p), n=length(y), case.weights = rep(1/n, n), reg_threshold, K=1000, p = dim(X)[2], start = NULL, constr_tol = 10^-4, sparse_warning=TRUE, maxit_SA = 1000, extra_step_NM=TRUE,
	short_cg_lim =15, long_cg_lim =500, ignore_below_zero=TRUE){


	if(any(is.null(dim(X)))) stop('X should be a matrix')
	if(is.null(start)) start <- rep(0,p+1)
	if(reg_threshold<=0) stop('invalid threshold')

	##### Set contant variables to have coefficient = 0, only analyze remaining variables in body of this function
		is_intercept <- apply(X,2,var)==0
		if(any(is_intercept) & sparse_warning) warning(sum(is_intercept), ' constant variables detected. Assigning these to have coefficient 0.')
		X_analyze <- X[,!is_intercept]
		reg_matrix_analyze <- reg_matrix[!is_intercept,!is_intercept]
		start_analyze <- start[c(TRUE,!is_intercept)]
		X_names <- colnames(X)

		rm(list=c('X','p','start'))
	#####

	invalid_constr <- TRUE 
	max_iter <- 20
	iter <- 1
	while(invalid_constr & iter < max_iter){
	
		fn_K <- function(w){
			get_hinge_obj(X=X_analyze,y=y,w=w, reg_threshold=reg_threshold, reg_matrix=reg_matrix_analyze, K=K, case.weights=case.weights)
		}
		gr_K <- function(w){
			get_hinge_grad(X=X_analyze,y=y,w=w, reg_threshold=reg_threshold, reg_matrix=reg_matrix_analyze, K=K, case.weights=case.weights)
		}
		
		# grad_based_step <- function(proposal, maxit_cg){Rcgmin(par=proposal, fn=fn_K, gr=gr_K,control=list(maxit=maxit_cg))}
		# grad_based_step <- function(proposal, maxit_cg){optim(par=proposal, fn=fn_K, gr=gr_K, method='BFGS',control=list(maxit=maxit_cg))}
		# grad_based_step <- function(proposal, maxit_cg){optim(par=proposal, fn=fn_K, gr=gr_K, method='CG',control=list(maxit=maxit_cg))}
		# grad_based_step <- function(proposal, maxit_cg){optimr(par=proposal, fn=fn_K, gr=gr_K, method='CG',control=list(maxit=maxit_cg))}
		grad_based_step <- function(proposal, maxit_cg){optimr(par=proposal, fn=fn_K, gr=gr_K, method='BFGS',control=list(maxit=maxit_cg))}
			#in ad hoc tests, this solver appeared to work best.
		

		soln <-
		warm_start <- grad_based_step(start_analyze, short_cg_lim)

		method <- 'grad-based'
		if(any(case.weights < 0)) method <- 'grad-based-plus-SA'

		if(method=='grad-based-plus-SA'){			

			parscale_SA <- pmax(
				max(abs(warm_start$par))*10^-8,
				abs(warm_start$par)
			)
			if(all(parscale_SA==parscale_SA[1])) parscale_SA <- rep(1,length(warm_start$par))
			### we would like SANN to give the same answer regardless of where it is run during a script -- change seed to ensure this.
			old_seed <- rnorm(1)*1000000
			SA_stop <- FALSE
			SA_fn <- function(proposal){
				#SA_stop is an ad hoc way to stop if we find a valid objective value below zero.
				if(SA_stop & ignore_below_zero) return(Inf)
				out <- fn_K(grad_based_step(proposal, short_cg_lim)$par)
				SA_stop <<- out < 0
				out
			}

			set.seed(0) 
				soln1_start_point <- optim(par=warm_start$par,
					fn=SA_fn, method='SANN',
					control=list(
						maxit=maxit_SA,
						parscale= parscale_SA))
			set.seed(old_seed)
			### 

			soln <-
			soln1 <- grad_based_step(proposal=soln1_start_point$par, maxit_cg=long_cg_lim)
		}
		if(method=='grad-based'){
			soln <-
			soln1 <- grad_based_step(proposal=start_analyze, maxit_cg=long_cg_lim)
		}

		

		if(extra_step_NM){
			nmk_error <- TRUE
			try({
				soln <- nmk(par=soln1$par, fn=fn_K)
				nmk_error <- FALSE
			})
			if(nmk_error){warning('nmk error, possibly due to instability.')}
		}

		out <- soln$par
		
		invalid_constr <- t(out[-1])%*%reg_matrix_analyze %*% out[-1] > reg_threshold * (1+constr_tol)
		K <- max(c(K,1000))*1000 #increase penalty if constraint is not met, and optimization must be re-run
		iter <- iter+1
	}
	if(invalid_constr) stop('Invalid constraint value')

	### Assign sparse categories to have 0 coefficient
	coef_vec <- c(NA,rep(0,length(is_intercept)))
	names(coef_vec) <- c('intercept', X_names)
	coef_vec[c(TRUE,!is_intercept)] <- out
	### 

	coef_vec
}





# Loop invariant statistics for hinge
get_suff_stats_linear_hinge <- function(reg_matrix, reg_threshold, start, ...){
	# compute loop invariant quantities for MCR binary search.
	fs <- get_full_sample(...)
	c(fs, list(reg_matrix=reg_matrix, start = start, reg_threshold=reg_threshold))
}

# regular loss for hinge
get_e0_hinge <- function(w,sst){
	mean(hinge_w(w,
		X=sst$full_X[!sst$is_perm,],
		y=sst$full_y[!sst$is_perm]))
}

# permuted loss for hinge
get_e1_hinge <- function(w,sst){
	mean(hinge_w(w,
		X=sst$full_X[sst$is_perm,],
		y=sst$full_y[sst$is_perm]))
}


