#Notes: 
# Throughout this code we refer to e_orig as e0, and to e_switch as e1.
# !! = important note to keep track of, or note for future changes.



#  _        __  __
# | |      |  \/  |
# | |      | \  / |
# | |      | |\/| |
# | |____  | |  | |
# |______| |_|  |_|




get_XtX_perm_lm <- function(X,p1){
	# Crossproduct of X-tilde, where X-tilde is the "full permuted" data frame.

	p <- dim(X)[2]
	n <- dim(X)[1]
	n1 <- rep(1,n)
	p2 <- setdiff(1:p,p1)
	XtX_perm <- crossprod(X) * (n-1)
	X1tX2_perm <- crossprod(X[,p1],n1) %*% crossprod(n1,(X[,p2])) - crossprod(X[,p1],X[,p2])
	XtX_perm[p1,p2] <- X1tX2_perm
	XtX_perm[p2,p1] <- t(X1tX2_perm)
	XtX_perm

	# Version that matches simplified Eq in paper, but is quadratic in n
	# p <- dim(X)[2]
	# n <- dim(X)[1]
	# p2_ind <- setdiff(1:p,p1)
	# W <- (1/(n-1)) * (tcrossprod(rep(1,n)) - diag(n)) #!! W written this way to match paper, although n-1 is cancelled out later
	# XtWX <- crossprod(X)
	# X1tWX2 <- t(X[,p1]) %*% W %*% X[,p2_ind]
	# XtWX[p1,p2_ind] <- X1tWX2
	# XtWX[p2_ind,p1] <- t(X1tWX2)
	# range(XtWX * (n-1) - XtX_perm)

}
get_ytX_perm_lm <- function(y,X,p1){
	# Crossproduct of y-tilde & X-tilde, which come from the "full permuted" dataset.
	# This is only meant to be used for linear models, not kernel models

	n <- length(y)
	n1 <- rep(1,n)
	ytX_perm <- crossprod(y,X)*(n-1)
	ytX_perm[1,p1] <- crossprod(y,n1) %*% crossprod(n1,(X[,p1])) - crossprod(y,X[,p1])

	ytX_perm 

	# Version that matches simplified Eq in paper, but is quadratic in n
	# n <- length(y)
	# W <- (1/(n-1)) * (tcrossprod(rep(1,n)) - diag(n)) #!! W written this way to match paper, although n-1 is cancelled out later
	# ytWX <- crossprod(y,X)
	# ytWX[1,p1] <- t(y) %*% W %*% X[,p1]
	# range(ytX_perm - ytWX * (n-1) )
}


# For internal use
# Get loop invariant statistics for MCR binary search with linear model class 
# and squared error loss
get_suff_stats_lm <- function(y,X,p1=NULL,tol=10^-10, reg_threshold=Inf, reg_matrix=NA){

	eigen_reg_matrix <- NA
	if(!any(is.na(reg_matrix))) eigen_reg_matrix <- eigen(reg_matrix)
	out <- list(
			tol=tol,
			n = length(y),
			yty_stnd = crossprod(y),
			ytX_stnd = crossprod(y,X),
			XtX_stnd = crossprod(X),
			reg_threshold=reg_threshold,
			reg_matrix=reg_matrix,
			eigen_reg_matrix = eigen_reg_matrix
		)
	if(!is.null(p1)){
		out$ytX_perm <- get_ytX_perm_lm(y=y, X=X, p1=p1)
		out$XtX_perm <- get_XtX_perm_lm(X=X, p1=p1)
	}
	out
}



##########
##########
##########
# Functions will also be used with suff_stats from get_suff_stats_kernel.

#' Get standard ("original") mean-squared-error for linear or kernel models. 
#' 
#' Can be used to calculate loss of a reference model, to determine epsilon.
#' @param model a coefficient vector (see \code{\link{fit_lm_regularized}})
#' @param suff_stats output from \code{\link{get_suff_stats_kernel}}, or from the internal function \code{get_suff_stats_lm}. See also \code{\link{precompute_mcr_objects_and_functions}}
#' @return Mean-squared-error in the sample.
#' @export
get_e0_lm <- function(model,suff_stats){
	# crossprod(y-X%*%model)/length(y)
	sst <- suff_stats
	out1 <- sst$yty_stnd - 2 * sst$ytX_stnd%*%model + t(model)%*% sst$XtX_stnd %*% model
	c(out1/sst$n)
}

# Get permuted loss for linear or kernel models
get_e1_lm <- function(model, suff_stats){
	sst <- suff_stats

	out1 <- sst$yty_stnd - 
		2 * sst$ytX_perm %*% model / (sst$n-1) + 
		t(model) %*% sst$XtX_perm %*% model / (sst$n-1)
	c(out1)/sst$n
}

# Get MR for linear or kernel models (internal)
get_MR_lm <- function(model, suff_stats){
	get_e1_lm(model=model, suff_stats)/get_e0_lm(model=model, suff_stats)
}

##########
##########
##########


#approximate workcheck
get_e1_lm_rand <- function(model, y, X, p1){ #only works for lm, not kernel.
	n <- length(y)
	X_rand <- X
	X_rand[,p1] <- X[sample(1:n),p1] # Approximate/biased, but close for large n
	crossprod(y-X_rand%*%model)/length(y)
}



# import MASS for ginv
#' @import MASS
get_h_min_lm_unregularized <- function(s, gam, loop_dep_args
=NULL, suff_stats){

	sst <- suff_stats

	if(!s %in% c(-1, 1)) stop('s input error')
	coef0 <- 1; coef1 <- 1
	if(s==-1) coef0 <- gam
	if(s==1) coef1 <- gam

	Q_gam <- coef0 * sst$XtX_stnd + coef1 * sst$XtX_perm / (sst$n-1)
	c_gam <- t(coef0 * sst$ytX_stnd + coef1 * sst$ytX_perm / (sst$n-1) )

	eigen_Q_gam <- eigen(Q_gam)
	if(any(eigen_Q_gam$values<0)) return(list(
		h_value = -Inf,
		e0_value  = -Inf,
		e1_value  = -Inf,
		beta_s_gam = NA
	))
	
	Q_gam_ginv <- ginv(Q_gam) # (!!) could do by hand to speed up
	beta_s_gam <- Q_gam_ginv %*% c_gam

	return(list(
		h_value = (1/sst$n) * c(
			sst$yty_stnd*(coef0+coef1) - t(c_gam) %*% Q_gam_ginv %*% c_gam
			), # "hat matrix" formula
		e0_value  = get_e0_lm(beta_s_gam,suff_stats=sst),
		e1_value  = get_e1_lm(beta_s_gam,suff_stats=sst),
		beta_s_gam = beta_s_gam
	))

}


#  _____    _       _
# |  __ \  (_)     | |
# | |__) |  _    __| |   __ _    ___
# |  _  /  | |  / _` |  / _` |  / _ \
# | | \ \  | | | (_| | | (_| | |  __/
# |_|  \_\ |_|  \__,_|  \__, |  \___|
#                        __/ |
#                       |___/

#' Can be used to fit either kernel models or regularized linear models
#' 
#' Minimizes \eqn{||Y_i - X \beta ||^2 +}\code{alpha}\eqn{\beta'(}\code{reg_matrix}\eqn{)\beta} across \eqn{\beta}. Alternatively, if \code{reg_threshold} is finite,  this function minimizes \eqn{||Y_i - X \beta ||^2} across \eqn{\beta} satisfying \eqn{\beta'(}\code{reg_matrix}\eqn{)\beta\le}\code{reg_threshold}.
#' 
#' @param suff_stats output from \code{\link{get_suff_stats_kernel}}, or from the internal function \code{get_suff_stats_lm}. See also \code{\link{precompute_mcr_objects_and_functions}}
#' @param alpha regularization parameter
#' @param reg_threshold regularization threshold
#' @param tol passed to \code{\link{solve_QP1QC}}
#' 
#' @export
#' @import qp1qc
#' @return The argmin over values of \eqn{\beta}.
fit_lm_regularized <- function(suff_stats, alpha=NA, reg_threshold=Inf, tol=NA){
	# ||Yi - XB ||^2 + alpha B'AB
	# y'y - 2y'X w + w'X'X w + alpha B'AB
	# deriv at
	# - 2 X'y + 2 X'X B_* + 2 alpha A B_* = 0
	# XtX B_* + alpha A B_* = X'y
	# w* = (XtX + alpha A)^-1 X'y

	#For kernel
	# ||Yi - sum_d k(X_i,D_d) w_d ||^2 + alpha w 'K_D w
	# ||Y  - K_stnd w ||^2 + alpha w'K_d w
	# y'y - 2y'K_stnd w + w'K_stnd'K_stnd w + alpha w'K_D w
	# deriv at
	# - 2 K_stnd'y + 2 K_stnd'K_stnd w* + 2 alpha K_D w = 0
	# K_stnd'K_stnd w* + alpha K_D w* = K_stnd'y
	# w* = (K_stnd'K_stnd + alpha K_D)^-1 K_stnd'y

	if(!mean(is.na(suff_stats)) %in% c(0,1)) warning('Partial NAs in suff_stats:',
		names(suff_stats)[which(is.na(suff_stats))])

	sst <- suff_stats
	if(is.infinite(reg_threshold)){
		if(is.na(alpha)) stop('alpha required when reg_threshold=Inf')
		B_star <- solve(sst$XtX_stnd + sst$reg_matrix*alpha*sst$n) %*% t(sst$ytX_stnd) #normalize alpha so that is not as affected when we increase sample size slightly from cross-validation to full sample. Otherwise, the size of XtX grows with n.
		return(B_star)
	}else{
		if(is.na(tol)) stop('tol must be set when reg_threshold<Inf')
		qp1qc_out <- solve_QP1QC(
			A_mat=sst$XtX_stnd, 
			a_vec = -2*c(sst$ytX_stnd), 
			B_mat =sst$reg_matrix, 
			b_vec = rep(0,length(sst$ytX_stnd)), 
			k=-reg_threshold, 
			tol=tol, verbose = FALSE)
		return(qp1qc_out$soln)
	}
}






#' Cross validate a kernel-least-squares, or kernel-regression model.
#' 
#' @export
#' @param y outcome vector
#' @param X covariate matrix
#' @param n_folds number of folds for CV
#' @param type either 'RKHS' or 'regression'
#' @param alpha penalization parameter for 'RKHS' (see \code{\link{fit_lm_regularized}})
#' @param dat_ref previously observed "reference" matrix of covariates, for 'RKHS'
#' @param warn_internal give warning if no dat_ref is not provided
#' @param kern_fun kernel function used
#' @param ... passed to \code{get_suff_stats_kernel}, for internal use.
CV_kernel <- function(y,X, n_folds=min(15,length(y)), type, alpha, dat_ref=NA, warn_internal=TRUE, kern_fun,...){
	#!! ... send to get_suff_stats, not fit_lm_regularized
	# (!!) Would be better if cv_mod allowed for fit_fun and report_fun to take indeces, optionally, instead of X,y.
		# Also would be worth considering trying to incorporate cv_mod here? (see also CV_lm)
		# This is currently not done to keep an eye on the de-meaning process,
		# but it should be possible to incorporate this also.
	# (!!) Can possibly improve speed by using X'X = sum_i X_i X_i'?
	
	if(any(is.na(dat_ref)) & warn_internal) warning('using internal data X as reference dataset (dat_ref)')
	
	get_suff_stat_fold <- function(ind, train, mu_y){ 
		if(any(is.na(dat_ref))){
			if( train) {dat_ref <- X[ind ,]}
			if(!train) {dat_ref <- X[!ind,]}
		}
		y_ind <- y[ind] - mu_y #!! demean
		X_ind <- X[ind,]
		if(sum(ind)==1) X_ind <- matrix(X_ind,nrow=1)
		get_suff_stats_kernel(y=y_ind, X=X_ind, reg_threshold=Inf, dat_ref=dat_ref, kern_fun=kern_fun, ...)
	}
	
	n_dat <- length(y)
	CV_err <- rep(NA,n_folds)
	CV_breaks <- cut(1:n_dat,breaks=n_folds)
	for(i in 1:n_folds){
		test_ind <- CV_breaks == levels(CV_breaks)[i]
		if(type=='RKHS'){
			mu_y_fold <- mean(y[!test_ind])
			sst_train_ind <- get_suff_stat_fold(ind=!test_ind, train = TRUE,  mu_y = mu_y_fold)
			sst_test_ind <-  get_suff_stat_fold(ind= test_ind, train = FALSE, mu_y = mu_y_fold)
			B_test_ind <- fit_lm_regularized(alpha=alpha, suff_stats=sst_train_ind)
			CV_err[i] <- get_e0_lm(suff_stats=sst_test_ind,model=B_test_ind)
		}
		if(type=='regression'){
			pred <- kernel_regression_prediction(
				X=X[test_ind,],
				X_ref=X[!test_ind,],
				y_ref=y[!test_ind], kern_fun=kern_fun)	
			CV_err[i] <- mean((y[test_ind]-pred)^2)
		}
	}
	mean(CV_err)
}




# !! Combine this with UNregularized function to reduce redundancy
get_h_min_lm_regularized <- function(s,gam,suff_stats,loop_dep_args
=NULL){


	sst <- suff_stats
	if(any(eigen(sst$reg_matrix)$value<=0)) stop('reg_matrix must be PD')


	if(!s %in% c(-1, 1)) stop('s input error')
	coef0 <- 1; coef1 <- 1
	if(s==-1) coef0 <- gam
	if(s==1) coef1 <- gam

	Q_gam <- coef0 * sst$XtX_stnd + coef1 * sst$XtX_perm / (sst$n-1)
	c_gam <- t(coef0 * sst$ytX_stnd + coef1 * sst$ytX_perm / (sst$n-1) )
	p<-length(sst$ytX_stnd)
	
	qp1qc_out <- solve_QP1QC(A_mat=Q_gam, a_vec = -2*c_gam, B_mat = sst$reg_matrix, b_vec = rep(0,p), k=-sst$reg_threshold, eigen_B_mat= sst$eigen_reg_matrix, tol=sst$tol, verbose=FALSE)
	beta_s_gam <- qp1qc_out$soln
	e0_value <- get_e0_lm(model=beta_s_gam,suff_stats=sst)
	e1_value <- get_e1_lm(model=beta_s_gam,suff_stats=sst)

	reg_value <- t(beta_s_gam)%*% sst$reg_matrix %*%beta_s_gam 
	if(reg_value > sst$reg_threshold * (1+sqrt( sst$tol))){
		warning(paste('Possible regularization error, magnitute ',signif(reg_value - sst$reg_threshold, digits=4)))
	}

	return(list(
		h_value = e0_value*coef0 + e1_value*coef1,
		e0_value  =e0_value,
		e1_value  =e1_value,
		reg_value = reg_value,
		beta_s_gam = beta_s_gam
	))

}



#  _                                   _
# | |                                 | |
# | | __   ___   _ __   _ __     ___  | |
# | |/ /  / _ \ | '__| | '_ \   / _ \ | |
# |   <  |  __/ | |    | | | | |  __/ | |
# |_|\_\  \___| |_|    |_| |_|  \___| |_|


#' Get loop-invariant quantities for empirical MCR binary search calculation, for regression in a RKHS. 
#' 
#' Creates `nrep_sample` copies of the sample, each of which permutes the covariate of interest in a different way. (Primarily intended for internal use.)
#' 
#' @param y scalar outcome vector
#' @param X matrix of covariates
#' @param kern_fun kernel function to use
#' @param dat_ref reference matrix of covariances
#' @param p1 index of X to permute
#' @param reg_threshold norm constraint for model class
#' @param tol tolerance used for calculations
#' @param verbose whether to report results during computation
#' @param nrep_sample see \code{\link{get_full_sample}}
#' @param warn_psd whether to warn that small numbers are added to make the kernel matrix positive definite
#' @param warn_duplicate whether to warn if duplicated rows of dat_ref are dropped
#' @param warn_dropped passed to \code{\link{get_full_sample}}
#' 
#' @importFrom kernlab kernelMatrix
#' @export
get_suff_stats_kernel <- function(y,X,dat_ref,kern_fun,p1=NULL, reg_threshold=Inf, tol=10^-10, nrep_sample=2, verbose=TRUE, warn_psd=TRUE, warn_duplicate=TRUE, warn_dropped=TRUE){
	# output will be identical to that of get_suff_stats_lm,
	# and can be used in exactly the same way, but is calculated differently (here)
	n<-dim(X)[1]
	if(any(duplicated(dat_ref))){
		if(warn_duplicate) warning('Dropping duplicate rows in dat_ref')
		dat_ref <- dat_ref[!duplicated(dat_ref),]
	}
	if(class(X)!='matrix' | class(dat_ref)!='matrix' ) stop('X and dat_ref must be of class "matrix"')
	if( #coefficient of variation test
		abs(mean(y))/
		sqrt( var(y)/length(y) + var(y)/dim(dat_ref)[1] )#sd of difference in iid sample means
		> 4) warning(paste('y may not be properly centered, check mean/intercept parameter used. Mean of y is ', signif(mean(y),3)))


	K_D <-   as.matrix(kernelMatrix(x=dat_ref,kernel=kern_fun))
	K_stnd <- as.matrix(kernelMatrix(x=X,y=dat_ref,kernel=kern_fun))

	eigen_K_D <- eigen(K_D)
	if(any(eigen_K_D$values<=tol)){
		#!! Perhaps there is a more elegant way to deal with this?
		err_K_D <- tail(eigen_K_D$values,3)
		add_to_K_D <- tol*2
		K_D <- K_D + diag(dim(K_D)[1]) * add_to_K_D
		if(warn_psd) warning(paste('singular K_D, last 3 eigenvalues are',
			paste(signif(err_K_D),collapse='; '),
			'. Adding',
			signif(add_to_K_D),'to diagonal'))
		eigen_K_D <- eigen(K_D)
	}
	if(any(eigen_K_D$values<=tol)) stop('Positive definite error')

	#Add intercept (not the same as in LM)
	col_vars <- apply(X, 2, var)
	if(any(col_vars==0)) stop('Constant column detected, which is not equivalent to adding an intercept. Consider dropping this column.')
	


	ytK_perm <- cp_K_perm <- NULL
	if(!is.null(p1)){
		
		###### 
		# Previous code, below, did not allow user-set nrep_sample, and did all permuations
		# which is not feasible for n > 5000 (as in our simulation population)
		#
			# The approach below of block calculations added together
			# is roughly 20% faster than storing the whole K_perm matrix
			#(for lenghth(y)=500; dim(dat_ref)=100)
			#
			# system.time({ 
			# 	cp_K_perm <- matrix(0,dim(dat_ref)[1],dim(dat_ref)[1])
			# 	ytK_perm <- matrix(0,1,dim(dat_ref)[1])
			
			# 	for(i in 1:n){
			# 		X_i <- X[-i,]
			# 		for(j in p1) X_i[,j] <- X[i,j] #!! technically not the same as in the manuscript [doing in one chunk can make issues with how it fills]
			# 		K_perm_i <- as.matrix(kernelMatrix(x=X_i,y=dat_ref,kernel=kern_fun))

			# 		cp_K_perm <- cp_K_perm + crossprod(K_perm_i)
			# 		ytK_perm <- ytK_perm + crossprod(y[-i],K_perm_i)
			# 	}
			# })
		######

		system.time({ 
			dim_cp <- dim(dat_ref)[1]
			cp_K_perm <- matrix(0,dim_cp,dim_cp)
			ytK_perm <- matrix(0,1,dim_cp)
			
			#!! In some cases, because of the observations left out of get_full_sample, the observations used here are different from those used in K_stnd (a warning indicates this event). A change does not seem necessary, but might be helpful.

			perm_ind <- get_full_sample(y=1:n,X=cbind(1:n,1:n),p1=1,n=n,nrep_sample=nrep_sample,seed=0, warn_dropped=warn_dropped)
			p1_perm_ind <- perm_ind$full_X[perm_ind$is_perm, 1]
			p2_perm_ind <- perm_ind$full_X[perm_ind$is_perm, 2]
			n_chunks <- nrep_sample-1
			n_perm_approx <- length(p2_perm_ind)
			if(n_chunks == 1) chunk_inds <- rep(1,n_perm_approx)
			if(n_chunks >  1) chunk_inds <- as.numeric(cut(1:n_perm_approx, breaks=n_chunks))

			if(verbose){
				cat('\n\n crossprod K_perm progress:\n')
				pb <- txtProgressBar(min=0,max=n_chunks)
			}
			for(i in 1:n_chunks){
				y_i <- y[p2_perm_ind[chunk_inds==i]]
				X_i <- X[p2_perm_ind[chunk_inds==i],]
				X_i[,p1] <- X[p1_perm_ind[chunk_inds==i],p1]
				K_perm_i <- as.matrix(kernelMatrix(x=X_i,y=dat_ref,kernel=kern_fun))


				cp_K_perm <- cp_K_perm + crossprod(K_perm_i) * (n*(n-1))/n_perm_approx
				ytK_perm  <- ytK_perm + crossprod(y_i,K_perm_i) * (n*(n-1))/n_perm_approx #n_perm_approx accounts for the fact that we are approximating a sume of n*(n-1)) terms with a sum of n_perm_approx terms.

				if(verbose) setTxtProgressBar(pb, i)

			}
		})
	}
	return(list(
		tol=tol,
		reg_threshold=reg_threshold,
		n=length(y),
		yty_stnd = crossprod(y),

		reg_matrix=K_D, 
		ytX_stnd =  crossprod(y,K_stnd),
		ytX_perm = ytK_perm,
		X_stnd = K_stnd,
		XtX_stnd = crossprod(K_stnd),
		XtX_perm = cp_K_perm,

		eigen_reg_matrix=eigen_K_D
	))
}

# Predict RKHS output
pred_RKHS <- function(K_stnd=NA, X=NA, dat_ref=NA,kern_fun=NA, model, ...){
	if(any(is.na(K_stnd))) K_stnd <- as.matrix(kernelMatrix(x=X,y=dat_ref,kernel=kern_fun))
	predictor <- K_stnd
	if( length(model) == dim(K_stnd)[2] + 1){
		predictor <- cbind(K_stnd,'intercept'=1) #add a column for the intercept
	}
	predictor %*% model
}

#' Get function norm for regression models in a RKHS
#' 
#' @param K_D kernel function applied to the reference data, in order to obtain a regularization matrix.
#' @param model model vector for which to calculate the norm.
#' @param dat_ref (optional) used to calculate \code{K_D}, if not specified.
#' @param kern_fun (optional) used to calculate \code{K_D}, if not specified.
#' @export
#' @return norm of \code{model} with respect to the regularization matrix \code{K_D}. This can be used in defining a model class, i.e., as the set of all models will norm less than that of a reference model.
norm_RKHS <- function(K_D=NA, dat_ref=NA, model, kern_fun){
	if(any(is.na(K_D))) K_D <-   as.matrix(kernelMatrix(x=dat_ref,kernel=kern_fun))
	model_no_intercept <- model[1:dim(K_D)[1]]
	if( abs(length(model) - dim(K_D)[1]) > 1 ) stop('dimension mismatch, even accounting for intercept.')
	t(model_no_intercept) %*% K_D %*% model_no_intercept
}

#' Get prediction from a (Nadaraya–Watson) kernel regression model
#' @param X covariate matrix for which predictions are desired
#' @param X_ref covariate matrix for existing data
#' @param y_ref corresponding output vector for existing data
#' @param kern_fun kernel function to use
kernel_regression_prediction <- function(X, X_ref, y_ref, kern_fun){
	#Simple kernel smoother estimator
	K <- as.matrix(kernelMatrix(x=X,y=X_ref, kernel=kern_fun))
	W <- K / tcrossprod(rowSums(K), rep(1,length(y_ref)))
	c(W %*% y_ref)

	# !! Note - See KRLS package to see how they handle binary covariates
		# However, this package does not appear to let you specify a present reference dataset; it always uses all your data.

}



