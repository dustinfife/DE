random.correlation <-
function(n){
	require(MASS)
	t <- mvrnorm(n,rep(0,n),diag(n))
	for (i in 1:n) {
		t[i,] <- t[i,]/sqrt(t(t[i,])%*%t[i,])
	}
	t%*%t(t)
}
