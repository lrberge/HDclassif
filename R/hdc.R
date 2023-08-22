
# TODO
# hdclassif_dim_choice => treat the case of cattell scree test with 2 EV

hdclassif_dim_choice <- function(ev, n, method, threshold, graph, noise.ctrl){
	# Selection of the intrinsic dimension 
	
	N <- sum(n)
	prop <- n/N
	K = ifelse(is.matrix(ev), nrow(ev), 1)
	
	# TODO
	# treat the case of cattell scree test with 2 EV
	
	if(graph){
		op = par(mfrow=c(K*(K<=3)+2*(K==4)+3*(K>4 && K<=9)+4*(K>9), 1+floor(K/4)-1*(K==12)+1*(K==7)))
		on.exit(par(op))
	}
	
	if(is.matrix(ev) && ncol(ev) <= 2){
		# d can be only equal to 1 if there are only 2 eigenvalues
		d = rep(1, K)
	} else if(is.matrix(ev) && K>1){
		p <- ncol(ev)
		if(method=="cattell"){
			dev <- abs(apply(ev, 1, diff))
			max_dev <- apply(dev, 2, max, na.rm=TRUE)
			dev <- dev/rep(max_dev, each=p-1)
			d <- apply((dev>threshold)*(1:(p-1))*t(ev[, -1]>noise.ctrl), 2, which.max)
			
			if(graph){
				
				for(i in 1:K){
					sub1 <- paste("Class #", i, ",  d", i, "=", d[i], sep="")
					
					Nmax <- max(which(ev[i, ]>noise.ctrl))-1
					plot(dev[1:(min(d[i]+5, Nmax)), i], type="l", col="blue", main=paste("Cattell's Scree-Test\n", sub1, sep=""), ylab=paste("threshold =", threshold), xlab="Dimension", ylim=c(0, 1.05))
					abline(h=threshold, lty=3) 	
					points(d[i], dev[d[i], i], col='red')
				}
				par(op)
			}
		} else if(method == "bic"){
			
			d <- rep(0, K)
			
			for (i in 1:K) {
				B <- c()
				
				Nmax <- max(max(which(ev[i, ]>noise.ctrl))-1, 1)
				p2 <- sum(!is.na(ev[i, ]))
				Bmax <- -Inf
				for (kdim in 1:Nmax){
					if ((d[i]!=0 & kdim>d[i]+10)) break
					a <- sum(ev[i, 1:kdim])/kdim
					b <- sum(ev[i, (kdim+1):p2])/(p2-kdim)
					if (b<0 | a<0){
						B[kdim] <- -Inf
					} else {
						L2 <- -1/2*(kdim*log(a)+(p2-kdim)*log(b)-2*log(prop[i])+p2*(1+1/2*log(2*pi))) * n[i]
						B[kdim] <- 2*L2 - (p2+kdim*(p2-(kdim+1)/2)+1) * log(n[i])
					}
					
					if ( B[kdim]>Bmax ){
						Bmax <- B[kdim]
						d[i] <- kdim
					}
				}
				
				if(graph){
					
					if(Nmax == 1){
						# we add some values
						ev_keep = ev[i, ][ev[i, ] > 0]
						p2 = length(ev_keep)
						for(kdim in 1:(min(10, length(ev_keep)))){
							a <- sum(ev[i, 1:kdim])/kdim
							b <- sum(ev[i, (kdim+1):p2])/(p2-kdim)
							if (b<0 | a<0){
								break
								B[kdim] <- -Inf
							} else {
								L2 <- -1/2*(kdim*log(a)+(p2-kdim)*log(b)-2*log(prop[i])+p2*(1+1/2*log(2*pi))) * n[i]
								B[kdim] <- 2*L2 - (p2+kdim*(p2-(kdim+1)/2)+1) * log(n[i])
							}
						}
						
					}
					
					plot(B, type='l', col=4, main=paste("class #", i, ",  d=", d[i], sep=''), ylab='BIC', xlab="Dimension")
					points(d[i], B[d[i]], col=2)
				}
			}
		}
	} else if(length(ev) <= 2){
		# idem, the number of intrinsic dimensions cannot be larger than 1 in the case of 2 eigenvalues
		d = 1
	} else {
		ev <- as.vector(ev)
		p <- length(ev)
		
		if(method == "cattell"){
			dvp <- abs(diff(ev))
			Nmax <- max(which(ev>noise.ctrl))-1
			if (p==2){
				d <- 1
			} else {
				d <- max(which(dvp[1:Nmax]>=threshold*max(dvp[1:Nmax])))
			} 
			
			diff_max <- max(dvp[1:Nmax])
			
			if(graph){
				plot(dvp[1:(min(d+5, p-1))]/diff_max, type="l", col="blue", main=paste("Cattell's Scree-Test\nd=", d, sep=''), ylab=paste("threshold =", threshold, sep=' '), xlab='Dimension', ylim=c(0, 1.05))
				abline(h=threshold, lty=3)	
				points(d, dvp[d]/diff_max, col='red')
			}
		} else if(method == "bic"){
			d <- 0
			Nmax <- max(max(which(ev>noise.ctrl))-1, 1)
			B <- c()
			Bmax <- -Inf
			for (kdim in 1:Nmax){
				
				if (d!=0 && kdim>d+10) break
				
				a <- sum(ev[1:kdim])/kdim
				b <- sum(ev[(kdim+1):p])/(p-kdim)
				if (b<=0 | a<=0){
					B[kdim] <- -Inf
				} else{
					L2 <- -1/2*(kdim*log(a)+(p-kdim)*log(b)+p*(1+1/2*log(2*pi)))*N
					B[kdim] <- 2*L2 - (p+kdim*(p-(kdim+1)/2)+1)*log(N)
				}
				
				if ( B[kdim]>Bmax ){
					Bmax <- B[kdim]
					d <- kdim
				}
			}
			
			if(graph){
				plot(B, type='l', col=4, main=paste("BIC criterion\nd=", d, sep=''), ylab='BIC', xlab="Dimension")
				points(d, B[d], col=2)
			}
		}
	}
	
	return(d)
}

hdclassif_bic <- function(par, p, data=NULL){
	model <- par$model
	K <- par$K
	d <- par$d
	b <- par$b
	a <- par$a
	mu <- par$mu
	N <- par$N
	prop <- par$prop
	
	if(length(b)==1){
		#update of b to set it as variable dimension models
		eps <- sum(prop*d)
		n_max <- if(model%in%c("ABQD", "AJBQD")) length(par$ev) else ncol(par$ev)
		b <- b*(n_max-eps)/(p-eps)
		b <- rep(b, length=K)
	}	
	
	if (length(a)==1){
		a <- matrix(a, K, max(d))
	} else if (length(a)==K) {
		a <- matrix(a, K, max(d))
	} else if (model=='AJBQD') {
		a <- matrix(a, K, d[1], byrow=TRUE)
	}
	
	if(min(a, na.rm=TRUE)<=0 | any(b<0)) return(-Inf)
	
	if (is.null(par$loglik)){
		som_a <- c()
		for (i in 1:K) som_a[i] <- sum(log(a[i, 1:d[i]]))
		L <- -1/2*sum(prop * (som_a + (p-d)*log(b) - 2*log(prop) + p*(1+log(2*pi))))*N
	} else if (model%in%c("ABQD", "AJBQD")){
		Q <- rep(list(par$Q), K)
		K_pen <- matrix(0, K, N)
		for (i in 1:K) {
			s <- sum(log(a[i, 1:d[i]]))
			X <- data-matrix(mu[i, ], N, p, byrow=TRUE)
			proj <- (X%*%Q[[i]])%*%t(Q[[i]])
			A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i, 1:d[i]], d[i]))
			B <- X-proj
			K_pen[i, ] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
		}
		A <- -1/2*t(K_pen)
		L <- sum(log(rowSums(exp(A-apply(A, 1, max))))+apply(A, 1, max))
	} else L <- par$loglik[length(par$loglik)]
	
	
	ro <- K*p+K-1
	tot <- sum(d*(p-(d+1)/2))
	D <- sum(d)
	d <- d[1]
	to <- d*(p-(d+1)/2)
	if (model=='AKJBKQKDK') m <- ro+tot+2*K+D
	else if (model=='AKBKQKDK') m <- ro+tot+3*K
	else if (model=='ABKQKDK') m <- ro+tot+2*K+1
	else if (model=='AKJBQKDK') m <- ro+tot+K+D+1
	else if (model=='AKBQKDK') m <- ro+tot+2*K+1
	else if (model=='ABQKDK') m <- ro+tot+K+2
	else if (model=='AKJBKQKD') m <- ro+K*(to+d+1)+1
	else if (model=='AKBKQKD') m <- ro+K*(to+2)+1
	else if (model=='ABKQKD') m <- ro+K*(to+1)+2
	else if (model=='AKJBQKD') m <- ro+K*(to+d)+2
	else if (model=='AKBQKD') m <- ro+K*(to+1)+2
	else if (model=='ABQKD') m <- ro+K*to+3
	else if (model=='AJBQD') m <- ro+to+d+2
	else if (model=='ABQD') m <- ro+to+3
	bic <- -(-2*L+m*log(N))
	
	#calcul ICL
	t = par$posterior
	if(!is.null(t)){
		# means we are in HDDC
		Z = ((t - apply(t, 1, max))==0) + 0
		icl = - (- bic - 2*sum(Z*log(t+1e-15)))
	} else {
		# Si HDDA, entropie est nulle => car classes pures
		icl = bic
	}
	
	return(list(bic = bic, icl = icl))
}

simuldata <- function(nlearn, ntest, p, K=3, prop=NULL, d=NULL, a=NULL, b=NULL){
	N=nlearn+ntest
	if (length(prop)==0) prop<-rep(1/K, K)
	else if (length(prop)!=K) stop("Proportions don't fit with the number of classes.")
	else prop<-prop/sum(prop)
	
	# Class sizes
	n<-floor(prop*N)
	N<-sum(n)
	
	#MEANS
	mu<-matrix(0, K, p)
	j<-sample(p, K)
	mu[cbind(1:K, j)]<-10
	
	# Intrinsic dimensions
	if ( length(d)==0 )	d<-sort(ceiling(runif(K, 0, 12*(p>20)+5*(p<=20 && p>=6)+(p<6)*(p-1))), decreasing=TRUE)
	else if ( length(d)!=K || !any(is.numeric(d)) ) stop("Wrong value of d.")
	
	# Orientation matrices
	Q<-vector(mode='list', length=K)
	for (i in 1:K) Q[[i]]<-qr.Q(qr(mvrnorm(p, mu=rep(0, p), Sigma=diag(1, p))))
	
	# Variance in the class-specific subspace
	if ( length(a)==0 ) a<-sort(ceiling(runif(K, 30, 350)))
	else if ( length(a)!=K || !any(is.numeric(a)) ) stop("Wrong value of a.")
	if ( length(b)==0 )b<-sort(ceiling(runif(K, 0, 25)))
	else if ( length(b)!=K || !any(is.numeric(b)) ) stop("Wrong value of b.")
	
	# Simulation
	S<-vector(mode='list', length=K)
	for (i in 1:K)	S[[i]]<-crossprod(Q[[i]]%*%sqrt(diag(c(rep(a[i], d[i]), rep(b[i], p-d[i])))))
	
	cls<-X<-NULL
	for (i in 1:K)	X<-rbind(X, mvrnorm(n[i], mu=mu[i, ], Sigma=S[[i]]))
	
	for (i in 1:K) cls<-c(cls, rep(i, n[i]))
	
	ind<-sample(1:N, N)
	prms<-list(a=a, b=b, prop=prop, d=d, mu=mu)
	data <- list(X=X[ind[1:nlearn], ], clx=cls[ind[1:nlearn]], Y=X[ind[(nlearn+1):N], ], cly=cls[ind[(nlearn+1):N]], prms=prms)
	
}

hdc_getComplexity = function(par, p){
	# Simple function to abtain the number of parameters of the models
	
	model <- par$model
	K <- par$K
	d <- par$d
	b <- par$b
	a <- par$a
	mu <- par$mu
	N <- par$N
	prop <- par$prop
	
	ro <- K*p+K-1
	tot <- sum(d*(p-(d+1)/2))
	D <- sum(d)
	d <- d[1]
	to <- d*(p-(d+1)/2)
	if (model=='AKJBKQKDK') m <- ro+tot+2*K+D
	else if (model=='AKBKQKDK') m <- ro+tot+3*K
	else if (model=='ABKQKDK') m <- ro+tot+2*K+1
	else if (model=='AKJBQKDK') m <- ro+tot+K+D+1
	else if (model=='AKBQKDK') m <- ro+tot+2*K+1
	else if (model=='ABQKDK') m <- ro+tot+K+2
	else if (model=='AKJBKQKD') m <- ro+K*(to+d+1)+1
	else if (model=='AKBKQKD') m <- ro+K*(to+2)+1
	else if (model=='ABKQKD') m <- ro+K*(to+1)+2
	else if (model=='AKJBQKD') m <- ro+K*(to+d)+2
	else if (model=='AKBQKD') m <- ro+K*(to+1)+2
	else if (model=='ABQKD') m <- ro+K*to+3
	else if (model=='AJBQD') m <- ro+to+d+2
	else if (model=='ABQD') m <- ro+to+3
	
	return(m)
}

hdc_getTheModel = function(model, all2models = FALSE){
	# Function used to get the models from number or names
	
	model_in = model
	
	if(!is.vector(model)) stop("The argument 'model' must be a vector.")
	
	if(anyNA(model)) stop("The argument 'model' must not contain any NA.")
	
	ModelNames <- c("AKJBKQKDK", "AKBKQKDK", "ABKQKDK", "AKJBQKDK", "AKBQKDK", "ABQKDK", "AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD", "AJBQD", "ABQD")
	
	model = toupper(model)
	
	if(length(model)==1 && model=="ALL"){
		if(all2models) model <- 1:14
		else return("ALL")
	}
	
	qui = which(model %in% 1:14)
	model[qui] = ModelNames[as.numeric(model[qui])]
	
	# We order the models properly
	qui = which(!model%in%ModelNames)
	if (length(qui)>0){
		if(length(qui)==1){
			msg = paste0("(e.g. ", model_in[qui], " is incorrect.)")	
		} else {
			msg = paste0("(e.g. ", paste0(model_in[qui[1:2]], collapse=", or "), " are incorrect.)")
		}
		stop("Invalid model name ", msg)
	}
	
	# warning:
	if(max(table(model))>1) warning("The model vector, argument 'model', is made unique (repeated values are not tolerated).")
	
	mod_num <- c()
	for(i in 1:length(model)) mod_num[i] <- which(model[i]==ModelNames)
	mod_num <- sort(unique(mod_num))
	model <- ModelNames[mod_num]
	
	return(model)
}

# We create a custom function to compute the eigen decomposition
hdc_myEigen = function(x, k, only.values=FALSE){
	
	if(k == ncol(x)){
		res = eigen(x, symmetric = TRUE, only.values = only.values)
	} else if (ncol(x) < 3){
		# The function eigs_sym can be used only for dataset with more than 3 vars
		res = eigen(x, symmetric = TRUE, only.values = only.values)
		# we put it in the right dimension
		res$values = res$values[1:k]
		if(!only.values){
			res$vector = res$vector[, 1:k]
		}
	} else {
		res = rARPACK::eigs_sym(x, k)
	}
	
	res
}

####
#### Utilities ####
####


addCommas = function(x) sapply(x, addCommas_single)

addCommas_single = function(x){
	# Cette fonction ajoute des virgules pour plus de
	# visibilite pour les (tres longues) valeurs de vraisemblance
	
	if(!is.finite(x)) return(as.character(x))
	
	s = sign(x)
	x = abs(x)
	
	decimal = x - floor(x)
	if(decimal>0) dec_string = substr(decimal, 2, 4)
	else dec_string = ""
	
	entier = as.character(floor(x))
	
	quoi = rev(strsplit(entier, "")[[1]])
	n = length(quoi)
	sol = c()
	for(i in 1:n){
		sol = c(sol, quoi[i])
		if(i%%3 == 0 && i!=n) sol = c(sol, ",")
	}
	
	res = paste0(ifelse(s==-1, "-", ""), paste0(rev(sol), collapse=""), dec_string)
	res
}

####
#### S3 methods ####
####


plot.hdc <- function(x, method=NULL, threshold=NULL, noise.ctrl=1e-8, ...){
	
	if(!is.null(method)){
		method = myAlerts(method, "method", "singleCharacterMatch.arg", "HDDC: ", c("cattell", "bic"))
	} else {
		if(!is.null(x$threshold) && is.numeric(threshold)){
			method = "cattell"
		} else {
			method = "bic"
		}
	}
	
	threshold <- if(!is.null(threshold)) threshold else if(!is.null(x$threshold) && is.numeric(threshold)) x$threshold else 0.2
	
	k <- x$K
	N <- x$N
	n <- x$prop*N
	
	if(is.null(x$com_ev)) d <- hdclassif_dim_choice(x$ev, n, method, threshold, TRUE, noise.ctrl)
	else d <- hdclassif_dim_choice(x$com_ev, n, method, threshold, TRUE, noise.ctrl)
}



#' Prediction method for \sQuote{hdc} class objects.
#' 
#' This function computes the class prediction of a dataset with respect to the model-based supervised and unsupervised classification methods \code{\link{hdda}} and \code{\link{hddc}}.
#' 
#' @method predict hdc
#'
#' @param object An \sQuote{hdc} class object obtained by using \code{\link{hdda}} or \code{\link{hddc}} function.
#' @param data A matrix or a data frame of observations, assuming the rows are the observations and the columns the variables. The data should be in the exact same format as the one that trained the model. Note that NAs are not allowed. 
#' @param cls A vector of the thue classes of each observation. It is optional and used to be compared to the predicted classes, default is NULL.
#' @param ... Not currently used.
#'
#' @return
#' \item{class}{vector of the predicted class.}
#' \item{prob}{The matrix of the probabilities to belong to a class for each observation and each class.}
#' \item{loglik}{The likelihood of the classification on the new data.}
#' If the initial class vector is given to the argument \sQuote{cls} then the adjusted rand index (ARI) is also returned. Also the following object is returned:
#' 	\item{ARI}{The confusion matrix of the classification.}
#'
#' @references 
#'Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High Dimensional Discriminant Analysis}, \emph{Communications in Statistics: Theory and Methods}, vol. \bold{36} (14), pp. 2607--2623
#'
#'Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High-Dimensional Data Clustering}, \emph{Computational Statistics and Data Analysis}, vol. \bold{52} (1), pp. 502--519
#'
#'Berge, L. Bouveyron, C. and Girard, S. (2012) \dQuote{HDclassif: An R Package for Model-Based 
#' Clustering and Discriminant Analysis of High-Dimensional Data}, \emph{Journal of Statistical Software}, 
#' \bold{46}(6), 1--29, url: \href{https://doi.org/10.18637/jss.v046.i06}{https://doi.org/10.18637/jss.v046.i06}
#'
#' @author
#'Laurent Berge, Charles Bouveyron and Stephane Girard
#'
#'
#' @seealso
#' The functions to do high dimensional classification \code{\link{hdda}} or clustering \code{\link{hddc}}.
#' 
#' @keywords hddc hdda clustering
#'
#' @examples
#'# Example 1:
#'data <- simuldata(1000, 1000, 50)
#'X <- data$X
#'clx <- data$clx
#'Y <- data$Y
#'cly <- data$cly
#'
#'#clustering of the gaussian dataset:
#'prms1 <- hddc(X, K=3, algo="CEM", init='param')      
#'           
#'#class vector obtained by the clustering:
#'prms1$class                   
#'
#'# only to see the good classification rate and 
#'# the Adjusted Rand Index:                     
#'res1 <- predict(prms1, X, clx)                                            
#'res2 <- predict(prms1, Y)       
#'
#'#the class predicted using hddc parameters on the test dataset:  
#'res2$class                                                           
#'
#'
#'# Example 2:
#'data(Crabs)
#'#clustering of the Crabs dataset:
#'prms3 <- hddc(Crabs[,-1], K=4, algo="EM", init='kmeans')        
#'res3 <- predict(prms3, Crabs[,-1], Crabs[,1])
#' 
#' 
#' 
#' 
predict.hdc <- function(object, data, cls=NULL, ...){
	#Extract variables:
	p <- ncol(data)
	N <- nrow(data)
	K <- object$K
	a <- object$a
	b <- object$b
	mu <- object$mu
	d <- object$d
	prop <- object$prop
	Q <- object$Q
	x <- as.matrix(data)
	confusion <- NULL
	ARI <- NULL
	
	if (length(N)==0) {
		N <- 1
		p <- length(data)
		x <- matrix(data, N, p)
	}
	
	if (length(object$scaling)!=0){
		x <- scale(x, center=object$scaling$mu, scale=object$scaling$sd)
	}
	
	if(length(b)==1){
		b <- rep(b, length=K)
	}
	
	if(length(a) == 1){
		a <- matrix(a, K, max(d))
	} else if(length(a) == K){
		a <- matrix(a, K, max(d))
	} else if(object$model == 'AJBQD'){
		a <- matrix(a, K, d[1], byrow=TRUE)
	}
	
	if(min(a, na.rm=TRUE)<=0 | min(b)<=0) stop("Some parameters A or B are negative. Prediction can't be done.\nThe reduction of the intrinsic dimensions or a more constrained model can be a solution.\nAlso, you can change the value of A's and B's manually by accessing the paramaters (though not recommended).\n", call.=FALSE)
	
	
	if(object$model=="AJBQD") {
		
		K_pen <- diag((mu%*%Q%*%diag(1/a[1, 1:d[1]], d[1]))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a[1, 1:d[1]], d[1]))%*%(t(Q)%*%t(x))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)+2*(mu%*%Q)%*%(t(Q)%*%t(x))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
		
	} else if (object$model=="ABQD") {
		
		K_pen <- diag(1/a[1]*(mu%*%Q)%*%(t(Q)%*%t(mu)))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))+2*(1/b[1]-1/a[1])*(mu%*%Q)%*%(t(Q)%*%t(x))
		
	} else {
		K_pen <- matrix(0, K, N)
		for (i in 1:K) {
			s <- sum(log(a[i, 1:d[i]]))
			X <- x-matrix(mu[i, ], N, p, byrow=TRUE)
			proj <- (X%*%Q[[i]])%*%t(Q[[i]])
			A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i, 1:d[i]], d[i]))
			B <- X-proj
			K_pen[i, ] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
		}
	}
	
	# The likelihood
	A <- -1/2*t(K_pen)
	A_max = apply(A, 1, max)
	loglik <- sum(log(rowSums(exp(A-A_max)))+A_max)
	
	t <- matrix(0, N, K, dimnames=list(1:N, 1:K))
	for (i in 1:K) t[, i] <- 1/rowSums(exp((K_pen[i, ]-t(K_pen))/2))
	result <- max.col(t)
	
	if (!is.null(object$kname)){
		result <- factor(result, labels=object$kname, levels=seq_along(object$kname))
		colnames(t) <- object$kname
	}
	
	if (!is.null(cls)){
		if(is.null(object$kname)){
			ARI <- hddc_ari(cls, result)
			cat("Adjusted Rand Index: ", ARI, ".\n", sep="")
		}
		else{
			confusion <- table(result, cls)
			dimnames(confusion) <- list('Predicted class'=object$kname, 'Initial class'=object$kname)
			cat("Correct classification rate: ", sum(diag(confusion))/N, ".\n", sep="")
			print(confusion)
		}
	}
	
	class(t) <- 'hd'
	
	list(class = result, posterior=t, loglik = loglik, confusion = confusion, ARI = ARI)
}

print.hd <- function(x, ...){
	class(x) <- NULL
	print.default(x, digits=3, na.print='.')
	class(x) <- 'hd'
}

print.hdc <- function(x, ...){
	if(length(x$kname)!=0) cat ("HIGH DIMENSIONAL DISCRIMINANT ANALYSIS\nMODEL: ", x$model, "\n", sep='')
	else cat ("HIGH DIMENSIONAL DATA CLUSTERING\nMODEL: ", x$model, "\n", sep='')
	print(x$prop)
	print(x$d)
	print(x$a)
	print(x$b)
	cat("BIC: ", x$BIC, "\n")
	if(!is.null(x$info)) cat(x$info, "\n")
	if(min(x$a, na.rm=TRUE)<0) cat("Information: a < 0\n")
	if(min(x$b)<10e-6) cat("Information: b < 10e-6\n")
}

####
#### PLOT ####
####

HDC_plot_criteria = function(res){
	
	all_crit = res$allCriteria
	all_crit = all_crit[order(all_crit$K, all_crit$model), ]
	K = all_crit$K
	BIC = all_crit$BIC
	model = all_crit$model
	model_unik = unique(model)
	nm = length(model_unik)
	
	if(length(unique(K)) == 1){
		plot(as.factor(model), BIC, ylab="BIC", main=paste("K =", K), xlab="model")
	} else {
		i = 0 
		for(m in model_unik){
			i = i + 1
			qui = which(model == m)
			quoi = all_crit[qui, ]
			if(m==model[1]){
				plot(quoi$K, quoi$BIC, ylab="BIC", ylim = range(all_crit$BIC, na.rm = TRUE))
			} else {
				lines(quoi$K, quoi$BIC, col = i, pch=i, lty=i, type="o", xlab="K")
			}
			legend(min(K, na.rm=TRUE), max(BIC), model_unik, col=1:nm, pch=1:nm, bty="n", lwd=1, cex=0.85, lty=1:nm)
		}
	}
	
}



####
#### Setters/Getters ####
####

#' Sets/gets the default 'show' argument in HDDC and HDDA
#'
#' Sets/gets the default value for 'show' argument in HDDC and HDDC. When \code{TRUE} then clustering information is returned at the end of the process. 
#'
#' @param show Single logical with default. Will specify the default value of the \code{show} argument in HDDA and HDDC.
#'
#'
#' @return
#' \code{getHDclassif.show} returns the default value.
#'
#' @examples
#'
#' data(Crabs)
#' 
#' # clustering of the Crabs dataset:
#' prms <- hddc(Crabs[,-1], K=4)  
#' # By default no information is displayed
#' 
#' # To show information:
#' prms <- hddc(Crabs[,-1], K=4, show = TRUE)  
#' 
#' # To set it permanently:
#' setHDclassif.show(TRUE)
#' prms <- hddc(Crabs[,-1], K=4)  
#' 
#' # to disable it permanently:
#' setHDclassif.show(FALSE)
#' 
#' 
#'
#'
setHDclassif.show = function(show){
	
	if(length(show) != 1 || !is.logical(show) || is.na(show)){
		stop("Argument 'show' must be a single logical.")
	}
	
	options("HDclassif.show" = show)
	
}

#' @rdname setHDclassif.show
"getHDclassif.show"

getHDclassif.show = function(){
	
	x = getOption("HDclassif.show")
	if(length(x) != 1 || !is.logical(x) || is.na(x)){
		stop("The value of getOption(\"HDclassif.show\") is currently not legal. Please use function setHDclassif.show to set it to an appropriate value.")
	}
	
	x
}



