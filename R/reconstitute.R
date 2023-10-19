reconstitute <- function(object){
	# fitted and residuals
	f <- stats::fitted(object) # involves centring and scaling up
	r <- stats::residuals(object, type="working")

	# Scale residuals if needed
	scal <- attr(r, "scaled:scale")

	# Scale and cetner
	if (!is.null(scal)) {
		r <- sweep(r, 2, scal, "*")
	}

	# Rebuild
	o <- f + sqrt(nrow(r)-1) * r

	# Return
	o
}
