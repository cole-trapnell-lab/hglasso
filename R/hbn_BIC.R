pseudo_ll <- function(x, theta) {
	n <- nrow(x)
	p <- ncol(X)
	part1 <- 0
	XtX <- t(x) %*% x
	for (j in 1:p) {
		for (j_prime in 1:p) {
			part1 <- part1 + theta[j, j_prime] * XtX[j, j_prime]
		}
	}

	part2 <- 0

	for (i in 1:n) {
		for (j in 1:p) {
			part2_b <- 0
			for (j_prime in 1:p) {
				if (j_prime != j) {
					part2_b <- part2_b + theta[j, j_prime] * x[i, j_prime]
				}
			}
			part2 <- part2 + log(1 + exp(theta[j, j] + part2_b))
		}
	}
	part1 - part2
}
hbnBIC <- function(hbn_object, theta, S = NULL) {
	n <- hbn_object$n	
	v <- length(hbn_object$hubind)	
	Zcard <- (sum(abs(hbn_object$Z)!=0)-hbn_object$p)/2
	Vcard <- hbn_object$V+t(x$V)
	diag(Vcard) <- 0
	Vcard <- sum(Vcard!=0)/2
	if (is.null(S)) S <- cov(x)

	return(list(BIC=-n*(pseudo_ll(x, theta))+sum(diag(S%*%x$Theta)))+log(n)*Zcard + log(n)*(v+c*(Vcard-v)))
}