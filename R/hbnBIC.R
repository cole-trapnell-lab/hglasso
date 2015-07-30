hbnBIC <- function(hbn_object, data, S = NULL, c=0.2) {
        n <- hbn_object$n
        v <- length(hbn_object$hubind)
        Zcard <- (sum(abs(hbn_object$Z)!=0)-hbn_object$p)/2
        Vcard <- hbn_object$V+t(hbn_object$V)
        diag(Vcard) <- 0
        Vcard <- sum(Vcard!=0)/2
        if (is.null(S)) S <- cov(data)
        return(list(BIC=-2*(pseudo_ll(data, hbn_object$Theta))+sum(diag(S%*%hbn_object$Theta))+log(n)*Zcard + log(n)*(v+c*(Vcard-v))))
}
