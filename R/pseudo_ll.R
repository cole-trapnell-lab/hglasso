pseudo_ll <- function(x, theta) {
	n <- nrow(x)
	p <- ncol(x)
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
