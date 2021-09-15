


demo_hddc  <- function(){
	myEnv <- new.env()
	data(list="Crabs",envir = myEnv)
	myData <- get("Crabs",envir=myEnv)
	devAskNewPage(ask = FALSE)
	X <- as.matrix(myData[,-1])
	clx <- myData[,1]
	algorithm <- c("EM","CEM","SEM")
	initialization <- c("kmeans","random","param")
	model <- c("AKBKQKDK","ABQKD","AJBQD")
	while (1){
		while(!any((algo <- readline("Choose the algorithm:\n1: EM; 2: CEM; 3: SEM; q or quit: exit \n"))==as.character(1:3))){
			if(any(tolower(algo)%in%c("q","quit"))) return(invisible())
		}
		
		while(!any((init <- readline("Choose initialization:\n1: kmeans; 2: random; 3: param; q or quit: exit \n"))==as.character(1:3))){
			if(any(tolower(init)%in%c("q","quit"))) return(invisible())
		}
		
		while(!any((mod <- readline("Choose the model:\n1: AkBkQkDk; 2: ABQkD; 3: AjBQD; q or quit: exit \n"))==as.character(1:3))){
			if(any(tolower(mod)%in%c("q","quit"))) return(invisible())
		}
		cat("hddc(data, classes, model=\"",model[as.numeric(mod)],"\",init=\"",initialization[as.numeric(init)],"\",algo=\"",algorithm[as.numeric(algo)],"\")\n",sep="")
		
		demo_hddc_crabs(X, 4, init=initialization[as.numeric(init)],algo=algorithm[as.numeric(algo)],model=model[as.numeric(mod)],ZZ=clx)
	}
}

demo_hddc_acp <- function(X, z, hd=NULL, ...){
	Q <- hd$Q
	MU <- hd$mu
	if (is.matrix(Q)) {
		svg <- vector(mode='list', length=4)
		for (i in 1:4) svg[[i]] <- as.matrix(Q)
		Q <- svg
	}
	
	X <- as.matrix(X)
	p <- ncol(X)
	k <- max(z)
	
	CO <- cov(X)
	Z <- -eigen(CO, symmetric=T)$vectors
	coul <- 1:4
	
	V <- X%*%Z
	patch <- c(3, 4, 8, 20)
	plot(V[, 1], V[, 2], col=z, pch=patch[z], ...)
	
	#drawing the projective space (a line); matrix(10, p, 1) is used only to have two points with the mean
	proj <- matrix(NA, k, p)
	for (i in 1:k) proj[i, ] <- tcrossprod(Q[[i]])%*%matrix(10, p, 1)+MU[i, ]
	
	x <- proj%*%Z
	y <- MU%*%Z
	points(y[, 1], y[, 2], col=coul[1:k], pch=19, lwd=7)
	
	for (i in 1:k) {
		pente <- (x[i, 2] - y[i, 2])/(x[i, 1] - y[i, 1])
		oo <- x[i, 2] - pente*x[i, 1]
		xb <- (2*y[i, 1] - sqrt(50^2/(pente^2+1)))/2
		xa <- (2*y[i, 1] + sqrt(50^2/(pente^2+1)))/2
		lines(c(xa, xb), oo+pente*c(xa, xb), col=coul[i], type='l')
	}
	
}

demo_hddc_crabs <- function(DATA, k=4, model='AKBKQKD',threshold=0.2, method='cattell',algo='EM',itermax=50, eps=1e-2, init='kmeans', ZZ=NULL, min.individuals=2, noise.ctrl=1e-8,...){ 
	com_dim <- 1
	Mod <- c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD")
	p <- ncol(DATA)
	N <- nrow(DATA)
	t <- matrix(0, N, k)
	if(model%in%Mod[7:14]) method <- 1
	if (init=='param'){
		MU <- colMeans(DATA)
		prop <- rep(1/k, k)
		S <- crossprod(DATA-matrix(MU, N, p, byrow=TRUE))/N
		donnees <- eigen(S, symmetric=TRUE)
		ev <- donnees$values
		d <- if(is.numeric(method)) method else hdclassif_dim_choice(ev, N, method, threshold, FALSE, noise.ctrl)
		a <- ev[1:d]
		b <- sum(ev[(d[1]+1):p])/(p-d[1])
		
		Q <- donnees$vectors[,1:d]
		mu <- mvrnorm(k, MU, S)
		
		K <- diag((mu%*%Q%*%diag(1/a, d, d))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a, d, d))%*%(t(Q)%*%t(DATA))+1/b*(diag(tcrossprod(mu))-2*mu%*%t(DATA)+2*(mu%*%Q)%*%(t(Q)%*%t(DATA))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
		
		t <- matrix(0, N, k)
		for (i in 1:k) t[,i]=1/rowSums(exp((K[i,]-t(K))/2))
	} else if (init=='kmeans') {
		mc <- match.call(expand.dots = FALSE)$...
		if (is.null(mc$algorithm)) alg="Hartigan-Wong"
		else alg=mc$algorithm
		if (is.null(mc$iter.max)) im=50
		else im=mc$iter.max
		if (is.null(mc$nstart)) nst=4
		else nst=mc$nstart
		cluster <- kmeans(DATA, k, iter.max=im, nstart=nst, algorithm=alg)$cluster
		for (i in 1:k) t[which(cluster==i),i] <- 1
	} else {
		t <- t(rmultinom(N, 1, rep(1/k, k)))
		compteur=1
		while(min(colSums(t))<1 && (compteur <- compteur+1)<5) t <- t(rmultinom(N, 1, rep(1/k, k)))
		if(min(colSums(t))<1) stop("Random initialization failed because of too many classes and too few observations")
	}
	
	likely <- c()
	I <- 0
	test <- Inf
	while ((I <- I+1)<=itermax && test>=eps){
		if (algo!='EM' && I!=1) t <- t2
		if (k>1 && (any(is.na(t)) || any(colSums(t>1/k)<min.individuals))) return(1)
		# browser()
		m <- hddc_m_step(DATA, k, t, model, threshold, method, noise.ctrl, 1, ncol(DATA))
		t <- hddc_e_step(DATA, m)
		
		L <- t$L
		t <- t$t
		if (algo=='CEM') {
			t2 <- matrix(0, N, k)
			t2[cbind(1:N, max.col(t))] <- 1
		}
		else if(algo=='SEM') { 
			t2 <- matrix(0, N, k)
			for (i in 1:N)	t2[i,] <- t(rmultinom(1, 1, t[i,]))
		}
		
		classes<-c()
		for (i in 1:N) classes[i]=which.max(t[i,])
		demo_hddc_acp(DATA, classes, m, xlab=paste('Iteration',I),ylab='',main="Clustering process",...)
		Sys.sleep(0.12)
		
		likely[I] <- L
		if (I!=1) test <- abs(likely[I]-likely[I-1])
	}
	
	cls <- max.col(t)
	ari = hddc_ari(ZZ, cls)
	cat("Adjusted Rand Index:",ari,"\n")
}



