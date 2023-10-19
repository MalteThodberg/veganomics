#### Variance decomposition ####

#' @describeIn plotDecomposition Plot total variance
#' @export
plotTotal <- function(object, color=NULL, shape=NULL, size=3){
	# Reconstitute if needed
	if(!is.null(object$CCA)){ #pCCA or CCA?
		message("Refitting PCA...")

		stopifnot(methods::is(object, "rda")) # not supported for CCA

		# Reconstitute
		# f <- fitted(object)
		# r <- residuals(object, type="working")
		# o <- f + sqrt(nrow(r)-1) * r

		# Rebuild input
		o <- reconstitute(object)

		# Redo PCA, scale if necessary
		if(is.null(attr(object$Ybar,"scaled:scale"))){
			object <- rda(o, scale=FALSE)
		}else{
			object <- rda(o, scale=TRUE)
		}

		# Set label
		header <- "Total: 100.0%"
	}else{
		header=NULL
	}

	# Scores
	o <- scores(object, display="sites", scaling = "sites")
	colnames(o) <- c("Dim1", "Dim2")
	o <- as.data.frame(o)

	# Axis labels
	tmp <- summary(eigenvals(object))
	dim1 <- paste0(colnames(tmp)[1], ": ", scales::percent(tmp[2,1], accuracy = 0.1))
	dim2 <- paste0(colnames(tmp)[2], ": ", scales::percent(tmp[2,2], accuracy = 0.1))

	# Plot
	ggplot(o, aes(x=Dim1, y=Dim2, color=color, shape=shape)) +
		geom_hline(yintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_vline(xintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_point(alpha=0.75, size=size) +
		labs(x=dim1,
				 y=dim2,
				 title=header) +
		coord_fixed()
}

#' @describeIn plotDecomposition Plot unconstrained variance
#' @export
plotUnconstrained <- function(object, color=NULL, shape=NULL, size=3){
	# Find rank
	r <- object$CCA$rank

	# Scores
	o <- scores(object, display=c("lc"), scaling = "sites", choices = c(r+1, r+2))
	colnames(o) <- c("Dim1", "Dim2")
	o <- as.data.frame(o)

	# Axis labels
	tmp <- summary(eigenvals(object, model = "unconstrained"))
	dim1 <- paste0(colnames(tmp)[1], ": ", scales::percent(tmp[2,1], accuracy = 0.1))
	dim2 <- paste0(colnames(tmp)[2], ": ", scales::percent(tmp[2,2], accuracy = 0.1))

	# Header
	antiR2 <- 1 - RsquareAdj(object)$r.squared
	header <- paste0("Unconstrained: ", scales::percent(antiR2, accuracy = 0.1))

	# Plot
	ggplot(o, aes(x=Dim1, y=Dim2, color=color, shape=shape)) +
		geom_hline(yintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_vline(xintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_point(alpha=0.75, size=size) +
		labs(x=dim1,
				 y=dim2,
				 title=header) +
		coord_fixed()
}

#' @describeIn plotDecomposition Plot constrained variance
#' @export
plotConstrained <- function(object, color=NULL, shape=NULL, size=3){
	# Scores
	o <- scores(object, display=c("lc", "wa"), scaling = "sites")

	# Find limies
	lims <- apply(rbind(o$sites, o$constraints), 2, range)

	# Format as data.frame
	o <- o$constraints
	colnames(o) <- c("Dim1", "Dim2")
	o <- as.data.frame(o)

	# Axis labels
	prefix <- switch(methods::is(object),
									 cca="CA",
									 rda="PC")
	tmp <- summary(eigenvals(object, model = "constrained"))
	dim1 <- paste0(prefix, ": ", scales::percent(tmp[2,1], accuracy = 0.1))

	if(ncol(tmp) > 1){
		dim2 <- paste0(prefix, ": ", scales::percent(tmp[2,2], accuracy = 0.1))
	}else{
		dim2 <- paste0(prefix, ": ~0%")
	}

	# Header
	R2 <- RsquareAdj(object)$r.squared
	header <- paste0("Constrained: ", scales::percent(R2, accuracy = 0.1))

	# Plot
	ggplot(o, aes(x=Dim1, y=Dim2, color=color, shape=shape)) +
		geom_hline(yintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_vline(xintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_point(alpha=0.75, size=size) +
		labs(x=dim1,
				 y=dim2,
				 title=header) +
		xlim(lims[1,1], lims[2,1]) +
		ylim(lims[1,2], lims[2,2]) +
		coord_fixed()
}

#' @describeIn plotDecomposition Plot RDA
#' @export
plotRDA <- function(object, color=NULL, shape=NULL, size=3){
	# Scores
	o <- scores(object, display=c("lc", "wa"), scaling = "sites")

	# Find limies
	lims <- apply(rbind(o$sites, o$constraints), 2, range)

	# Format as data.frame
	o <- o$sites
	colnames(o) <- c("Dim1", "Dim2")
	o <- as.data.frame(o)

	# Header
	R2 <- RsquareAdj(object)$adj.r.squared
	header <- paste0("Adjusted: ", scales::percent(R2, accuracy = 0.1))

	# Plot
	ggplot(o, aes(x=Dim1, y=Dim2, color=color, shape=shape)) +
		geom_hline(yintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_vline(xintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_point(alpha=0.75, size=size) +
		labs(x="RDA1",
				 y="RDA2",
				 title=header) +
		xlim(lims[1,1], lims[2,1]) +
		ylim(lims[1,2], lims[2,2]) +
		coord_fixed()
}


#' Visualize RDA or CCA with variance partitioning
#'
#' @param object rda or cca object
#' @param color vector
#' @param shape vector
#' @param size scalar
#' @param color_scale Scale object from ggplot2
#' @param shape_scale Scale object from ggplot2
#'
#' @return ggplot
#' @export
plotDecomposition <- function(object, color=NULL, shape=NULL, size=3, color_scale=NULL, shape_scale=NULL){
	# Unconstrained, constrained, total and RDA
	plot_list <- list(total=plotTotal(object, color=color, shape=shape),
										rda=plotRDA(object, color=color, shape=shape),
										constrained=plotConstrained(object, color=color, shape=shape),
										unconstrained=plotUnconstrained(object, color=color, shape=shape))

	# Add color scale
	if(!is.null(color_scale)){
		stopifnot(methods::is(color_scale, "Scale"))
		plot_list <- lapply(plot_list, function(i) i + color_scale)
	}

	# Add shape scale
	if(!is.null(shape_scale)){
		stopifnot(methods::is(shape_scale, "Scale"))
		plot_list <- lapply(plot_list, function(i) i + shape_scale)
	}

	# Meta header
	if(is.null(object$pCCA)){
		header <- "Redundancy Analysis (RDA)"
	}else{
		header <- "Partial Redundancy Analysis (pRDA)"
	}

	wrap_plots(plot_list, ncol=2, nrow=2) +
		plot_layout(guides="collect", widths = 1) +
		plot_annotation(title = header)
}

#### Introspection plots ####

#' Visualize RDA or CCA with study design
#' @inheritParams plotDecomposition
#' @param bp_col character: Color or biplot colors
#' @param bp_scale scalar: scaling factor of biplot scale
#' @export
plotDesign <- function(object, color=NULL, shape=NULL, size=3, bp_col="black", bp_scale=0.9){
	# Scores
	o <- scores(object, display=c("wa", "bp"), scaling = "sites")
	bp <- as.data.frame(o$biplot)
	o <- as.data.frame(o$sites)

	# Biplot scaling
	mult <- min(
		(max(o[,2]) - min(o[,2])/(max(bp[,2])-min(bp[,2]))),
		(max(o[,1]) - min(o[,1])/(max(bp[,1])-min(bp[,1]))))

	bp <- bp*mult

	# Header
	R2 <- RsquareAdj(object)$adj.r.squared
	header <- paste0("Adjusted: ", scales::percent(R2, accuracy = 0.1))

	# Plot
	ggplot(o) +
		geom_hline(yintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_vline(xintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_point(aes(x=RDA1, y=RDA2, color=color, shape=shape), alpha=0.75, size=size) +
		geom_point(data=bp, aes(x=RDA1, y=RDA2), size=size, color=bp_col) +
		geom_segment(data=bp, aes(x=0, xend=RDA1,
															y=0, yend=RDA2),
								 color=bp_col) +
		ggrepel::geom_label_repel(data=bp, aes(x=RDA1, y=RDA2, label=rownames(bp)),
															color=bp_col) +
		labs(x="RDA1",
				 y="RDA2",
				 title=header) +
		coord_fixed()
}

#' Visualize RDA or CCA design strain (spider plot)
#' @inheritParams plotDecomposition
#' @export
plotStrain <- function(object, color=NULL, shape=NULL, size=3){
	# Extract scores
	o <- scores(object, display=c("lc", "wa"), scaling = "sites")

	# Format as to-from
	o <- cbind(o$sites, o$constraints)
	o <- as.data.frame(o)
	colnames(o) <- c("sites1", "sites2", "constraints1", "constraints2")

	# Explained variance
	R2 <- RsquareAdj(object)$adj.r.squared
	header <- paste0("Adjusted: ", scales::percent(R2, accuracy = 0.1))

	# Plot
	ggplot(o, aes(x=sites1, y=sites2, xend=constraints1, yend=constraints2, color=color)) +
		geom_hline(yintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_vline(xintercept = 0, linetype="dashed", color="black", alpha=0.25) +
		geom_segment(alpha=0.5) +
		geom_point(aes(shape=shape), alpha=0.75, size=size) +
		labs(x="RDA1", y="RDA2", title=header) +
		coord_fixed()
}

#### Heap ####

# plotDesign <- function(object, color=NULL, shape=NULL, size=3, bp_col="black", bp_scale=0.9, ...){
# 	# Scores
# 	o <- scores(object, display=c("wa", "bp"), scaling = "sites")
# 	bp <- as.data.frame(o$biplot * 5)
# 	o <- as.data.frame(o$sites)
#
# 	# Scores
# 	bp$RDA1 <- bp$RDA1 * (max(o$RDA1) / max(bp$RDA1))
# 	bp$RDA2 <- bp$RDA2 * (max(o$RDA2) / max(bp$RDA2))
#
# 	# Header
# 	R2 <- RsquareAdj(object)$adj.r.squared
# 	header <- paste0("Adjusted: ", scales::percent(R2, accuracy = 0.1))
#
# 	# Plot
# 	ggplot(o) +
# 		geom_hline(yintercept = 0, linetype="dashed", color="black", alpha=0.25) +
# 		geom_vline(xintercept = 0, linetype="dashed", color="black", alpha=0.25) +
# 		geom_text(data=bp, aes(x=RDA1, y=RDA2, label=rownames(bp)),
# 							alpha=0.50,
# 							color=bp_col) +
# 		geom_segment(data=bp, aes(x=0, xend=RDA1 * bp_scale,
# 															y=0, yend=RDA2 * bp_scale),
# 								 alpha=0.50,
# 								 color=bp_col,
# 								 arrow = arrow(...)) +
# 		geom_point(aes(x=RDA1, y=RDA2, color=color, shape=shape), alpha=0.75, size=size) +
# 		labs(x="RDA1",
# 				 y="RDA2",
# 				 title=header) +
# 		coord_fixed()
# }



