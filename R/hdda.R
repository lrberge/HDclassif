



#' High Dimensional Discriminant Analysis
#' 
#' HDDA is a model-based discriminant analysis method assuming each class of the dataset live in a proper Gaussian subspace which is much smaller than the original one, the hdda.learn function calculates the parameters of each subspace in order to predict  the class of new observation of this kind.
#'
#' @inheritParams hddc
#' @param cls The vector of the class of each observations, its type can be numeric or string.
#' @param d_select Either \dQuote{Cattell} (default), \dQuote{BIC} or \dQuote{CV}. See details for more information. This parameter selects which method to use to select the intrinsic dimensions.
#' @param graph It is for comparison sake only, when several estimations are run at the same time (either when using several models, or when using cross-validation to select the best dimension/threshold). If graph = TRUE, the plot of the results of all estimations is displayed. Default is FALSE.
#' @param cv.dim A vector of integers. Only when d_select=\dQuote{CV}. Gives the dimensions for which the CV is to be done. Note that if some dimensions are greater than what it is possible to have, those are taken off.
#' @param cv.threshold A vector of floats strictly within 0 and 1. Only when d_select=\dQuote{CV}. Gives the thresholds for which the CV is to be done.
#' @param cv.vfold An integer. Only when d_select=\dQuote{CV}. It gives the number of different subsamples in which the dataset is split. If \dQuote{cv.vfold} is greater than the number of observations, then the program equalize them.
#' @param LOO If TRUE, it returns the results (classes and posterior probabilities) for leave-one-out cross-validation.
#' 
#' @details
#' Some information on the signification of the model names:
#' \describe{
#' 	\item{Akj are the parameters of the classes subspaces:}{
#' 		\itemize{
#' 		\item{if Akj: each class has its parameters and there is one parameter for each dimension}
#' 		\item{if Ak: the classes have different parameters but there is only one per class}
#' 		\item{if Aj: all the classes have the same parameters for each dimension (it's a particular case with a common orientation matrix)}
#' 		\item{if A: all classes have the same one parameter}
#' 		}
#' 	}
#' 
#' 	\item{Bk are the noises of the classes subspaces:}{
#' 		\itemize{
#' 			\item{If Bk: each class has its proper noise}
#' 			\item{if B:  all classes have the same noise}
#' 		}
#' 	}
#' 
#' 	\item{Qk is the orientation matrix of each class:}{ 
#' 		\itemize{
#' 			\item{if Qk: all classes have its proper orientation matrix}
#' 			\item{if Q: all classes have the same orientation matrix}
#' 		}
#' 	}
#' 	
#' 	\item{Dk is the intrinsic dimension of each class:}{ 
#' 		\itemize{
#' 			\item{if Dk: the dimensions are free and proper to each class}
#' 			\item{if D: the dimension is common to all classes}
#' 		}
#' 	}
#' }
#' The model \dQuote{all} will compute all the models, give their BIC and keep the model with the highest BIC value.
#' Instead of writing the model names, they can also be specified using an integer.  1 represents the most general model (\dQuote{AkjBkQkDk}) while 14 is the most constrained (\dQuote{ABQD}), the others  number/name matching are given below. Note also that several models can be run at once, by using a vector of models (e.g. model = c("AKBKQKD","AKJBQKDK","AJBQD") is equivalent to model = c(8,4,13); to run the 6 first models, use model=1:6). If all the models are to be run, model="all" is faster than model=1:14. 
#' \tabular{lcclc}{
#' AkjBkQkDk \tab   1   \tab   \tab  AkjBkQkD \tab   7   \cr 
#' AkBkQkDk \tab   2   \tab \tab  AkBkQkD \tab   8   \cr   
#' ABkQkDk \tab   3   \tab  \tab ABkQkD \tab   9   \cr   
#' AkjBQkDk \tab   4   \tab  \tab  AkjBQkD \tab   10   \cr   
#' AkBQkDk \tab   5   \tab  \tab  AkBQkD \tab   11   \cr   
#' ABQkDk \tab   6   \tab  \tab  ABQkD \tab   12  \cr
#' AjBQD \tab 13 \tab  \tab ABQD \tab 14
#' }
#' 
#' The parameter d, is used to select the intrinsic dimensions of the subclasses. Here are his definictions:
#' 		\itemize{
#' 			\item{\dQuote{Cattell}:}{
#' 				The Cattell's scree-test is used to gather the intrinsic dimension of each class. If the model is of common dimension (models 7 to 14), the scree-test is done on the covariance matrix of the whole dataset.
#' 			}
#' 			\item{\dQuote{BIC}:}{
#' 				The intrinsic dimensions are selected with the BIC criterion. See Bouveyron \emph{et al.} (2010) for a discussion of this topic.
#' 				For common dimension models, the procedure is done on the covariance matrix of the whole dataset.
#' 			}
#' 			\item{\dQuote{CV}:}{
#' 				A V-fold cross-validation (CV) can be done in order to select the best threshold (for all models) or the best common dimensions (models 7 to 14).  
#' 				The V-fold cross-validation is done for each dimension (respectively threshold) in the argument \dQuote{cv.dim} (resp. \dQuote{cv.threshold}), then the dimension (resp. threshold) that gives the best good classification rate is kept.  
#' 				The dataset is split in \dQuote{cv.vfold} (default is 10) \emph{random} subsamples, then CV is done for each sample: each of them is used as validation data while the remaining data is used as training data. For sure, if \dQuote{cv.vfold} equals the number of observations, then this CV is equivalent to a leave-one-out.
#' 			}
#' 		}
#'
#' @return
#' hdda returns an 'hdc' object; it's a list containing:
#' \item{ model }{The name of the model.}
#' \item{ k }{The number of classes.}
#' \item{ d }{The dimensions of each class.}
#' \item{ a }{The parameters of each class subspace.}
#' \item{ b }{The noise of each class subspace.}
#' \item{ mu }{The mean of each variable for each class.}
#' \item{ prop }{The proportion of each class.}
#' \item{ ev }{The eigen values of the var/covar matrix.}
#' \item{ Q }{The orthogonal matrix of orientation of each class.}
#' \item{ kname }{The name of each class.}
#' \item{ BIC }{The BIC value of the model used.}
#' \item{ scaling }{The centers and the standard deviation of the original dataset.}
#' 
#' @references
#' Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High Dimensional Discriminant Analysis}, \emph{Communications in Statistics: Theory and Methods}, vol. \bold{36} (14), pp. 2607--2623
#' 
#' Bouveyron, C. Celeux, G. and Girard, S. (2010) \dQuote{Intrinsic dimension estimation by maximum likelihood in probabilistic PCA}, Technical Report 440372, Universite Paris 1 Pantheon-Sorbonne
#' 
#' Berge, L. Bouveyron, C. and Girard, S. (2012) \dQuote{HDclassif: An R Package for Model-Based Clustering and Discriminant Analysis of High-Dimensional Data}, \emph{Journal of Statistical Software}, \bold{46}(6), 1--29, url: \href{http://www.jstatsoft.org/v46/i06/}{http://www.jstatsoft.org/v46/i06/}
#' 
#' @author
#' Laurent Berge, Charles Bouveyron and Stephane Girard
#' 
#' @seealso
#' \code{\link{hddc}}, \code{\link{predict.hdc}}, \code{\link{plot.hdc}}
#'
#' @examples
#' # Example 1:
#' data<-simuldata(1000, 1000, 50, K=5)
#' X <- data$X
#' clx <- data$clx
#' Y <- data$Y
#' cly <- data$cly
#' # we get the HDDA parameters:
#' prms1 <- hdda(X, clx)         
#' 
#' cl1 <- predict(prms1, Y, cly)
#' # the class vector of Y estimated with HDDA:
#' cl1$class                     
#' 
#' # another model is used:
#' prms1 <- hdda(X, clx, model=12)
#' #model=12 is equivalent to model="ABQkD"     
#' cl1 <- predict(prms1, Y, cly) 
#' 
#' # Example 2:
#' data(wine)
#' a <- wine[,-1]
#' z <- wine[,1]
#' prms2 <- hdda(a, z, model='all', scaling=TRUE, d_select="bic", graph=TRUE)
#' cl2 <- predict(prms2, a, z)
#' 
#' # getting the best dimension
#' # using a common dimension model
#' # we do LOO-CV using cv.vfold=nrow(a)
#' prms3 <- hdda(a, z, model="akjbkqkd", d_select="CV", cv.vfold=nrow(a), scaling=TRUE, graph=TRUE)
#' 
#' cl3 <- predict(prms3, a, z)
#' 
#' # Example 3:
#' # Validation with LOO
#' prms4 = hdda(a, z, LOO=TRUE, scaling=TRUE)
#' sum(prms4$class==z) / length(z)
#' 
hdda <- function(data, cls, model='AkjBkQkDk', graph=FALSE, d_select="Cattell", threshold=0.2, com_dim=NULL, show=getHDclassif.show(), scaling=FALSE, cv.dim=1:10, cv.threshold=c(.001, .005, .05, 1:9*0.1), cv.vfold=10, LOO=FALSE, noise.ctrl=1e-8, d){
	
	# For compatibility with old versions of HDclassif
	if(!missing(d) & missing(d_select)) d_select = d
	
	#
	# CONTROLS
	#
	
	call = match.call()
	hdda_control(call)
	# Control of match.args:
	d_select = myAlerts(d_select, "d_select", "singleCharacterMatch.arg", "HDDA: ", c("cattell", "bic", "cv"))
	# We get the model names, properly ordered
	model = hdc_getTheModel(model)
	
	ModelNames <- c("AKJBKQKDK", "AKBKQKDK", "ABKQKDK", "AKJBQKDK", "AKBQKDK", "ABQKDK", "AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD", "AJBQD", "ABQD", "ALL")
	Mod2 <- c("AKJBKQKDK", "AKBKQKDK ", "ABKQKDK ", "AKJBQKDK ", "AKBQKDK ", "ABQKDK  ", "AKJBKQKD ", "AKBKQKD ", "ABKQKD  ", "AKJBQKD ", "AKBQKD  ", "ABQKD  ", "AJBQD  ", "ABQD   ")
	
	cls <- as.factor(cls)
	names <- levels(cls)
	z <- unclass(cls)
	data <- as.matrix(data)
	K <- max(z)
	
	if (scaling) {
		data <- scale(data)
		scaling <- list(mu=attr(data, "scaled:center"), sd=attr(data, "scaled:scale"))
	} else scaling <- NULL
	
	N <- nrow(data)
	p <- ncol(data)
	n <- table(z)
	
	#
	# Leave one out
	#
	
	if(LOO){
		
		if(model%in%ModelNames[1:14] && !is.null(com_dim)) d_select <- com_dim
		
		class <- rep(NA, N)
		posterior <- matrix(NA, N, K)
		for(i in 1:N){
			prms <- NULL
			try(prms <- hdda_prms(data[-i, ], z[-i], model, threshold, d_select, names, noise.ctrl), silent=TRUE)
			if(!is.null(prms)) {
				res <- NULL
				try(res <- predict.hdc(prms, data[i, ]), silent=TRUE)
				if(!is.null(res)){
					class[i] <- res$class
					posterior[i, ] <- res$posterior
				}
			}
		}
		class <- factor(class, labels=names, levels=seq_along(names))
		return(list(class=class, posterior=posterior))
	}
	
	
	#
	# Cross Validation
	#
	
	if(d_select=="cv"){
		
		d_select = "cattell"
		
		d.max <- if(model%in%ModelNames[7:12]) min(n, p)-1 else min(N, p)-1
		
		cv.dim <- sort(cv.dim, decreasing=TRUE)
		cv.dim <- cv.dim[cv.dim<=d.max]
		if(length(cv.dim)==0) stop("cv.dim must be an integer stricly inferior to the dimension.", call.=FALSE)
		
		cv.threshold <- sort(cv.threshold) 
		cv.threshold <- cv.threshold[cv.threshold>=0 & cv.threshold<=1]
		if(length(cv.threshold)==0) stop("cv.threshold must be a float within 0 and 1.\n", call.=FALSE)
		cv.vfold <- if(cv.vfold<N) cv.vfold else N
		
		u <- sample(1:N)
		ind <- c()
		for(i in 1:cv.vfold) ind[i] <- if(i==1) floor(N/cv.vfold) else floor((N-sum(ind))/(cv.vfold+1-i))
		fin <- cumsum(ind)
		debut <- c(1, fin[-cv.vfold]+1)
		
		if(model%in%ModelNames[7:14]) {
			n_cv <- length(cv.dim)
			cv.threshold <- rep(.5, n_cv)
		} else{
			n_cv <- length(cv.threshold)
			cv.dim <- rep("cattell", n_cv)
		}
		
		res <- fails <- rep(0, n_cv)
		N2 <- rep(N, n_cv)
		for(j in 1:cv.vfold){
			ind <- u[debut[j]:fin[j]]
			prms <- NULL
			i <- 0
			while((i <- i+1)<=n_cv && is.null(prms) ){
				
				try(prms <- hdda_prms(data[-ind, ], z[-ind], model, cv.threshold[i], cv.dim[i], names, noise.ctrl), silent=TRUE)
				
				if(!is.null(prms)){
					try(res[i] <- res[i] + sum(predict.hdc(prms, data[ind, ])$class==cls[ind]), silent=TRUE)
				} else {
					N2[i] <- N2[i]-length(ind)
					fails[i] <- fails[i] + 1
				}
			}
			
			if(i<=n_cv) for(i in i:n_cv){
				
				if(model%in%ModelNames[1:6]){
					d <- hdclassif_dim_choice(prms$ev, as.vector(table(z[-ind])), "cattell", cv.threshold[i], FALSE, noise.ctrl)
				} else d <- rep(cv.dim[i], K)
				
				
				if(model%in%ModelNames[13:14]){
					prms$Q <- prms$Q[, 1:d[1]]
				} else {
					for (ii in 1:K) {
						if (prms$d[ii] > 1) {
							prms$Q[[ii]] <- prms$Q[[ii]][, 1:d[ii]]
						}
					}
				} 
				
				prms$d <- d	
				prms_bis <- hdda_prms_bis(model, prms, p)
				try(res[i] <- res[i] + sum(predict.hdc(prms_bis, data[ind, ])$class==cls[ind]), silent=TRUE)
			}
		}
		
		if(show){
			if(model%in%ModelNames[7:14]) cat("\t Model  \t dim\t CV\n")
			else cat("\t Model  \tthreshold\t CV\n")
			for(i in n_cv:1){
				if(model%in%ModelNames[7:14]) cat('\t', Mod2[model==ModelNames], '\t', cv.dim[i], "\t", res[i]/N2[i]*100, if(fails[i]>0) paste(" Info: failed", fails, "times"), '\n')
				else cat('\t', Mod2[model==ModelNames], '\t', cv.threshold[i], "\t\t", res[i]/N2[i]*100, if(fails[i]>0) paste(" Info: failed", fails, "times"), '\n')
			}
		}
		
		res <- res/N2*100
		res <- res[n_cv:1]
		cv.dim <- cv.dim[n_cv:1]
		cv.threshold <- cv.threshold[n_cv:1]
		
		if(model%in%ModelNames[7:14]){
			d <- com_dim <- cv.dim[which.max(res)]
			if(show) cat("Best dimension with respect to the CV results: ", d, ".\n", sep="")
			if(graph){
				barplot(res-100/K, names.arg=cv.dim, offset=100/K, col="blue", xlab="Dimensions", ylab="Correct classification rate", axes=FALSE, main=paste("Cross-Validation\n(chosen dim=", d, ")", sep=""))
				axis(2, at=floor(100/K+(max(res)-100/K)/5*0:5))
			}
			d_select = d # a number 
		} else{
			d_select <- "cattell"
			threshold <- cv.threshold[which.max(res)]
			if(show) cat("Best threshold with respect to the CV results: ", threshold, ".\n", sep="")
			if(graph){
				barplot(res-100/K, names.arg=cv.threshold, offset=100/K, col="blue", xlab="Thresholds", ylab="Correct classification rate", axes=FALSE, main=paste("Cross-Validation\nthreshold=", threshold, sep=""))
				axis(2, at=floor(100/K+(max(res)-100/K)/5*0:5))
			}
		}
	}	
	
	#
	# Any model
	#
	
	if(length(model)>1){
		nm = length(model)
		e <- vector(mode="list", length=nm)
		BIC <- ICL <- c()
		
		for(i in 1:nm){
			e[[i]] <- hdda_prms(data, z, model[i], threshold, d_select, names, noise.ctrl, com_dim)
			BIC[i] <- hdclassif_bic(e[[i]], p)$bic
			ICL[i] <- hdclassif_bic(e[[i]], p)$icl
		}
		
		prms <- e[[which.max(BIC)]]
		prms$BIC <- max(BIC, na.rm=TRUE)
		prms$scaling <- scaling
		prms$threshold <- threshold
		
		if(show){
			cat(" # :\t Model \t   BIC\n")
			for(i in 1:nm){
				if(i<10) cat(' ')
				wng <- if(any(e[[i]]$b<10e-6) | any(e[[i]]$a<10e-6, na.rm=TRUE)) "info: b < 10e-6" else ""
				cat(i, ':\t', Mod2[ModelNames == model[i]], '\t', BIC[i], wng, '\n')
			}
			cat("\nSELECTED: Model ", prms$model, ", BIC=", prms$BIC, ".\n", sep="")
		}
		
		if(graph){			
			BIC <- BIC[!is.na(BIC)]
			min_b=min(BIC[BIC!=-Inf])
			max_b=max(BIC, na.rm=TRUE)
			BIC[BIC==-Inf] <- min_b
			barplot(BIC-min_b, names.arg=model, offset=min_b, col="blue", xlab="models", ylab="BIC", axes=FALSE, main=paste("BIC for all models\n(chosen model=", prms$model, ")", sep=""))
			axis(2, at=floor(min_b+(max_b-min_b)/5*0:5))
		}
		
		class(prms) <- 'hdc'
		return(prms)
	} else if(model=="ALL"){
		e <- vector(mode="list", length=14)
		BIC <- ICL <- c()
		
		#models with var dim
		e[[1]] <- hdda_prms(data, z, ModelNames[1], threshold, d_select, names, noise.ctrl)
		for (i in 2:6) e[[i]] <- hdda_prms_bis(ModelNames[i], e[[1]], p)
		
		#models with common dim	
		e[[7]] <- hdda_prms(data, z, ModelNames[7], threshold, d_select, names, noise.ctrl, com_dim)
		for (i in 8:12) e[[i]] <- hdda_prms_bis(ModelNames[i], e[[7]], p)
		
		#models 13 and 14: common var/covar matrix
		e[[13]] <- hdda_prms(data, z, ModelNames[13], threshold, d_select, names, noise.ctrl, com_dim)
		e[[14]] <- hdda_prms_bis(ModelNames[14], e[[13]], p)
		
		#BIC calculation
		for(i in 1:14){
			BIC[i] <- hdclassif_bic(e[[i]], p)$bic
			ICL[i] <- hdclassif_bic(e[[i]], p)$icl
		}
		
		prms <- e[[which.max(BIC)]]
		prms$BIC <- max(BIC, na.rm=TRUE)
		prms$scaling <- scaling
		prms$threshold <- threshold
		
		if(show){
			cat(" # :\t Model \t   BIC\n")
			for(i in 1:14){
				if(i<10) cat(' ')
				wng <- if(any(e[[i]]$b<10e-6) | any(e[[i]]$a<10e-6, na.rm=TRUE)) "info: b < 10e-6" else ""
				cat(i, ':\t', Mod2[i], '\t', BIC[i], wng, '\n')
			}
			cat("\nSELECTED: Model ", prms$model, ", BIC=", prms$BIC, ".\n", sep="")
		}
		
		if(graph){
			min_b <- min(BIC[BIC!=-Inf])
			max_b <- max(BIC)
			BIC[BIC==-Inf] <- min_b
			barplot(BIC-min_b, names.arg=1:14, offset=min_b, col="blue", xlab="models", ylab="BIC", axes=FALSE, main=paste("BIC for all models\n(chosen model=", prms$model, ")", sep=""))
			axis(2, at=floor(min_b+(max_b-min_b)/5*0:5))
		}
		
		class(prms) <- 'hdc'
		return(prms)
	} else {
		prms <- hdda_prms(data, z, model, threshold, d_select, names, noise.ctrl, com_dim)
		prms$BIC <- hdclassif_bic(prms, p)$bic
		prms$ICL <- hdclassif_bic(prms, p)$icl
		prms$scaling <- scaling
		prms$threshold <- threshold
		
		class(prms) <- 'hdc'
		return(prms)
	}
}

hdda_prms <- function(data, cls, model, threshold, method, kname, noise.ctrl, com_dim=NULL){
	
	# Some parameters
	p <- ncol(data)
	N <- nrow(data)
	K <- max(cls)
	com_ev <- NULL
	info <- NULL
	n <- as.vector(table(cls))
	prop <- matrix(n/N, 1, K, dimnames=list(c(''), "Prior probabilities of groups:"=kname))
	
	mu <- matrix(rowsum(data, cls)/n, K, p, dimnames=list("Class"=kname, "Group means:"=paste('V', 1:p, sep='')))
	
	#Calculation of Var/covar matrices and of eigenvectors
	
	if( model%in%c("AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD", "AJBQD", "ABQD") ){
		# For common dimension models => we first get the common eigenvalues
		
		if (N<p) {
			Y <- matrix(0, N, p)
			for (i in 1:K){
				qui = which(cls == i)
				Y[qui, ] <- (data[qui, ] - matrix(mu[i, ], length(qui), p, byrow=TRUE))/sqrt(N)
			}
			
			if(model%in%c("AJBQD", "ABQD")){
				donnees <- eigen(tcrossprod(Y), symmetric=TRUE)
			} else {
				donnees <- eigen(tcrossprod(Y), symmetric=TRUE, only.values=TRUE)
			}
		} else {
			W <- matrix(0, p, p)
			for (i in 1:K){
				W <- W + prop[i]*crossprod(data[which(cls==i), ] - matrix(mu[i, ], sum(cls==i), p, byrow=TRUE))/n[i]
			}
			
			if(model%in%c("AJBQD", "ABQD")){
				donnees <- eigen(W, symmetric=TRUE)
			} else {
				donnees <- eigen(W, symmetric=TRUE, only.values=TRUE)
			} 
		}	
		ev <- com_ev <- donnees$values
	}
	
	if(!model%in%c("AJBQD", "ABQD")){
		if(any(n<p)){
			Y <- vector(mode='list', length=K)
		}
		
		Q <- vector(mode='list', length=K)
		
		ev <- matrix(NA, K, min(max(n), p))
		
		for(i in which(n<p)){
			Y[[i]] <- (data[which(cls==i), ] - matrix(mu[i, ], sum(cls==i), p, byrow=TRUE))/sqrt(n[i])
			donnees <- eigen(tcrossprod(Y[[i]]), symmetric=TRUE)
			ev[i, 1:n[i]] <- donnees$values
			Q[[i]] <- donnees$vectors
		}
		
		for(i in which(n>=p)){
			donnees <- eigen(crossprod(data[which(cls==i), ] - matrix(mu[i, ], sum(cls==i), p, byrow=TRUE))/n[i], symmetric=TRUE)
			ev[i, ] <- donnees$values
			Q[[i]] <- donnees$vectors
		}
	}
	
	#Intrinsic dimension calculations + graphical display
	
	if(model%in%c("AJBQD", "ABQD")){
		
		if(!is.null(com_dim)) method <- com_dim
		
		if(method%in%c("cattell", "bic")) method <- hdclassif_dim_choice(com_ev, n, method, threshold, FALSE, noise.ctrl)
		d <- rep(method, K)
		
	} else if (model%in%c('AKJBKQKD', 'AKBKQKD', 'ABKQKD', 'AKJBQKD', 'AKBQKD', 'ABQKD')){
		
		if(!is.null(com_dim)) method <- com_dim
		
		if(method%in%c("cattell", "bic")) method <- hdclassif_dim_choice(com_ev, n, method, threshold, FALSE, noise.ctrl)
		
		d <- rep(method, K)
		if( d[1]>min(n, p)-1 ) {
			d[] <- min(n, p)-1
			info <- paste("Information: d has been lowered to", d[1], "because of the class", kname[which.min(n)], "which has", min(n), "observations.")
		}
		dmax <- if(any(ev<noise.ctrl, na.rm=TRUE)) max(min(unlist(apply(ev<noise.ctrl, 1, which)))-2, 1) else Inf
		if(d[1] > dmax) d[] <- dmax
	} else {
		# browser()
		d <- hdclassif_dim_choice(ev, n, method, threshold, FALSE, noise.ctrl)
	} 
	
	# Setup of Qi matrices
	# Dans le cas ou N<p, on normalise les matrices
	
	if (model%in%c("AJBQD", "ABQD")){
		if (N>=p) {
			Q <- matrix(donnees$vectors[, 1:d[1]], p, d[1])
		} else {
			Q <- matrix(t(Y)%*%donnees$vectors[, 1:d[1]], p, d[1])
			normalise <- c()
			for(i in 1:d[1]) normalise[i] <- as.double(crossprod(Q[, i]))
			Q <- Q/matrix(sqrt(normalise), p, d[1], byrow=TRUE)
		}
	} else{
		for(i in which(n>=p)){
			Q[[i]] <- matrix(Q[[i]][, 1:d[i]], p, d[i])
		}
		
		for(i in which(n<p)){
			Q[[i]] <- t(Y[[i]])%*%(Q[[i]][, 1:d[i]])
			normalise <- c()
			for (j in 1:d[i]) normalise[j] <- as.double(crossprod(as.matrix(Q[[i]][, j])))
			Q[[i]] <- Q[[i]]/matrix(sqrt(normalise), p, d[i], byrow=TRUE)
		}
	}
	
	# Calculation of the remaining parameters of the selected model
	
	# PARAMETER a
	if ( model%in%c('AKJBKQKDK', 'AKJBQKDK', 'AKJBKQKD', 'AKJBQKD') ){
		ai <- matrix(NA, K, max(d), dimnames=list("Class"=kname, "Akj:"=paste("a", 1:max(d), sep='')))
		for (i in 1:K) ai[i, 1:d[i]] <- ev[i, 1:d[i]]
	} else if ( model%in%c('AKBKQKDK', 'AKBQKDK' , 'AKBKQKD', 'AKBQKD') ){
		ai <- matrix(NA, 1, K, dimnames=list(c("Ak:"), kname))
		for (i in 1:K) ai[i] <- sum(ev[i, 1:d[i]])/d[i]
	} else if (model=="AJBQD"){
		ai <- matrix(ev[1:d[1]], 1, d[1], dimnames=list(c("Aj:"), paste('a', 1:d[1], sep='')))
	} else if (model=="ABQD"){
		ai <- matrix(sum(ev[1:d[1]])/d[1], dimnames=list(c("A:"), c('')))
	} else {
		a <- 0
		eps <- sum(prop*d)
		for (i in 1:K) a <- a + sum(ev[i, 1:d[i]])*prop[i]
		ai <- matrix(a/eps, dimnames=list(c("A:"), c('')))
	}
	
	# PARAMETER b
	if ( model%in%c('AKJBKQKDK', 'AKBKQKDK', 'ABKQKDK', 'AKJBKQKD', 'AKBKQKD', 'ABKQKD') ){
		bi <- matrix(NA, 1, K, dimnames=list(c("Bk:"), kname))
		for(i in which(n>=p)) bi[i] <- sum(ev[i, (d[i]+1):p])/(p-d[i])
		for(i in which(n<p)) bi[i] <- sum(ev[i, (d[i]+1):n[i]])/(p-d[i])
	} else if ( model%in%c("ABQD", "AJBQD") ){
		if (N>=p) bi <- matrix(sum(ev[(d[1]+1):p])/(p-d[1]), dimnames=list(c("B:"), c('')))
		else bi <- matrix(sum(ev[(d[1]+1):N])/(N-d[1]), dimnames=list(c("B:"), c('')))
	} else{
		b <- 0
		eps <- sum(prop*d)
		for(i in which(n>=p)) b <- b + sum(ev[i, (d[i]+1):p])*prop[i]
		for(i in which(n<p)) b <- b + sum(ev[i, (d[i]+1):n[i]])*prop[i]
		bi <- matrix(b/(min(max(n), p)-eps), dimnames=list(c("B:"), c('')))
	}
	
	d <- matrix(d, 1, K, dimnames=list(c('dim:'), "Intrinsic dimensions of the classes:"=kname))
	
	class(prop) <- class(mu) <- class(ai) <- class(bi) <- class(d) <- class(ev) <- "hd"
	
	list(model=model, K=K, d=d, a=ai, b=bi, mu=mu, prop=prop, ev=ev, Q=Q, kname=kname, info=info, N=N, com_ev=com_ev)
}

hdda_prms_bis <- function(model, par, p){
	# Used to update the parameters a and b when we already have the eigenvectors and the Qs
	
	N <- par$N
	K <- par$K
	ev <- par$ev
	d <- par$d
	kname <- par$kname
	prop <- par$prop
	n <- prop*N
	
	# PARAMETER a
	if ( model%in%c('AKJBKQKDK', 'AKJBQKDK', 'AKJBKQKD', 'AKJBQKD') ){
		ai <- matrix(NA, K, max(d), dimnames=list("Class"=kname, "Akj:"=paste("a", 1:max(d), sep='')))
		for (i in 1:K) ai[i, 1:d[i]] <- ev[i, 1:d[i]]
	} else if ( model%in%c('AKBKQKDK', 'AKBQKDK' , 'AKBKQKD', 'AKBQKD') ){
		ai <- matrix(NA, 1, K, dimnames=list(c("Ak:"), kname))
		for (i in 1:K) ai[i] <- sum(ev[i, 1:d[i]])/d[i]
	} else if (model=="AJBQD"){
		ai <- matrix(ev[1:d[1]], 1, d[1], dimnames=list(c("Aj:"), paste('a', 1:d[1], sep='')))
	} else if (model=="ABQD"){
		ai <- matrix(sum(ev[1:d[1]])/d[1], dimnames=list(c("A:"), c('')))
	} else {
		a <- 0
		eps <- sum(prop*d)
		for (i in 1:K) a <- a + sum(ev[i, 1:d[i]])*prop[i]
		ai <- matrix(a/eps, dimnames=list(c("A:"), c('')))
	}
	
	# PARAMETER b
	if ( model%in%c('AKJBKQKDK', 'AKBKQKDK', 'ABKQKDK', 'AKJBKQKD', 'AKBKQKD', 'ABKQKD') ){
		bi <- matrix(NA, 1, K, dimnames=list(c("Bk:"), kname))
		for(i in which(n>=p)) bi[i] <- sum(ev[i, (d[i]+1):p])/(p-d[i])
		for(i in which(n<p)) bi[i] <- sum(ev[i, (d[i]+1):n[i]])/(p-d[i])
	} else if ( model%in%c("ABQD", "AJBQD") ){
		if (N>=p) bi <- matrix(sum(ev[(d[1]+1):p])/(p-d[1]), dimnames=list(c("B:"), c('')))
		else bi <- matrix(sum(ev[(d[1]+1):N])/(N-d[1]), dimnames=list(c("B:"), c('')))
	} else{
		b <- 0
		eps <- sum(prop*d)
		for(i in which(n>=p)) b <- b + sum(ev[i, (d[i]+1):p])*prop[i]
		for(i in which(n<p)) b <- b + sum(ev[i, (d[i]+1):n[i]])*prop[i]
		bi <- matrix(b/(min(max(n), p)-eps), dimnames=list(c("B:"), c('')))
	}
	
	class(ai) <- class(bi) <- "hd"
	
	list(model=model, K=K, d=d, a=ai, b=bi, mu=par$mu, prop=par$prop, ev=ev, Q=par$Q, kname=par$kname, info=par$info, N=N, com_ev=par$com_ev)
}

####
#### CONTROLS ####
####

hdda_control = function(call){
	
	prefix = "HDDA: "
	myCallAlerts(call, "data", "matrix, data.frame", 3, TRUE, prefix)
	myCallAlerts(call, "cls", "(integer, factor)Vector", 3, TRUE, prefix)
	myCallAlerts(call, "model", "vector", 3, FALSE, prefix)
	myCallAlerts(call, "graph", "singleLogical", 3, FALSE, prefix)
	myCallAlerts(call, "d_select", "singleCharacter", 3, FALSE, prefix)
	myCallAlerts(call, "threshold", "numericVectorGT0LT1", 3, FALSE, prefix)
	myCallAlerts(call, "com_dim", "singleIntegerGE1", 3, FALSE, prefix)
	myCallAlerts(call, "show", "singleLogical", 3, FALSE, prefix)
	myCallAlerts(call, "scaling", "singleLogical", 3, FALSE, prefix)
	myCallAlerts(call, "cv.dim", "integerVectorGE1", 3, FALSE, prefix)
	myCallAlerts(call, "cv.threshold", "numericVectorGE0LE1", 3, FALSE, prefix)
	myCallAlerts(call, "cv.vfold", "singleIntegerGE1", 3, FALSE, prefix)
	myCallAlerts(call, "LOO", "singleLogical", 3, FALSE, prefix)
	myCallAlerts(call, "noise.ctrl", "singleNumericGE0", 3, FALSE, prefix)
	
	myCallAlerts(call, "mc.cores", "singleIntegerGE1", 3, FALSE, prefix)
	
	####
	#### SPECIFIC controls 
	####
	
	# Getting some elements
	data = eval.parent(call[["data"]], 2)
	cls = eval.parent(call[["cls"]], 2)
	
	# No NA in the data:
	if (any(is.na(data))) stop("NA values in the data are not supported. Please remove them beforehand.")
	if (any(is.na(cls))) stop("NA values in the classes (argument 'cls') are not supported. Please remove them beforehand.")
	
	# length must match
	if (nrow(data)!=length(cls)) stop ("The class vector does not fit the data (they are of different lengths)", call.=FALSE)
	
	# model controls:
	model = eval.parent(call[["model"]], 2)
	d_select = eval.parent(call[["d_select"]], 2)
	if(!is.null(model) & !is.null(d_select)){
		if(d_select == "CV" && (length(model)>1 || (length(model)==1 && toupper(model)=="ALL"))) stop("A specific model must be chosen for choosing the threshold/dimension with cross-validation.\n")
	}
	
	# Common dimension models
	com_dim = eval.parent(call[["com_dim"]], 2)
	cls <- as.factor(cls)
	names <- levels(cls)
	n = as.vector(table(cls))
	N = sum(n)
	p = ncol(data)
	ModelNames <- c("AKJBKQKDK", "AKBKQKDK", "ABKQKDK", "AKJBQKDK", "AKBQKDK", "ABQKDK", "AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD", "AJBQD", "ABQD", "ALL")
	if( (any(model%in%ModelNames[7:12]) && !is.null(com_dim) && com_dim > min(n, p)-1) ) stop("com_dim has to be lower or equal to ", min(n, p)-1, if(p>min(n)) paste(" because of the class", names[which.min(n)], "which has", min(n), "observations.\n") else paste(" because there are only", p, "dimensions.\n"))
	if(any(model%in%ModelNames[13:15]) && !is.null(com_dim) && com_dim > min(N, p)-1) stop("com_dim has to be lower or equal to ", min(N, p)-1, if(p<=N) paste(" because there are only", p, "dimensions.\n") else paste(" because there are only", N, "observations.\n"))
	
	# LOO
	LOO = eval.parent(call[["LOO"]], 2)
	if(!is.null(LOO) && LOO){
		if(!is.null(d_select) && d_select=="CV" && is.null(com_dim)) stop("The cross validation method to select the dimensions cannot be used when doing LOO. Otherwise, you can select manually the dimension for common dimension models using 'com_dim'.\n")
		if(!is.null(model) && (length(model)>1 || model=="ALL")) stop("A specific model must be chosen for LOO.\n")
	}
	
}
















