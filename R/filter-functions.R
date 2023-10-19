#' Filter matrix for low counts
#'
#' Wrapper for model.matrix and filterByExpr for removing low counts genes from expression matrix
#'
#' @inheritParams stats::model.matrix
#' @inheritParams edgeR::filterByExpr
#' @return filter and transposed count matrix
#' @export
filterCounts <- function(object, data,
												 y,
												 min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7){

	# Setup up
	mod <- stats::model.matrix(object, data=data)

	# Filter
	above <- edgeR::filterByExpr(y=y,
												design=mod,
												min.count = min.count,
												min.total.count = min.total.count,
												large.n = large.n,
												min.prop = min.prop)
	o <- y[above,]

	# Return
	t(o)
}

#' Normalize and log-transform count matrix
#'
#' Wrapper for calcNormFactors and cpm from the edgeR package.
#' @inheritParams edgeR::DGEList
#' @inheritParams edgeR::cpm
#' @param ... additional arguments passed to edgeR::calcNormFactors
#'
#' @return transposed and log-normalized expression matrix
#' @export
normalizeCounts <- function(counts, log=TRUE, prior.count=5, ...){
	# Transpose back to RNA
	o <- t(counts)

	# Normalize
	o <- edgeR::DGEList(o)
	o <- edgeR::calcNormFactors(o, ...)
	o <- edgeR::cpm(o, log=TRUE, prior.count = prior.count)

	# Return
	t(o)
}

#' Variance-stabilize count matrix
#'
#' Wrapper for varianceStabilizingTransformation from the DESeq2 package.
#'
#' @inheritParams DESeq2::DESeqDataSetFromMatrix
#' @inheritParams DESeq2::varianceStabilizingTransformation
#'
#' @return transposed and log-normalized expression matrix
#' @export
stabilizeCounts <- function(design, colData, countData, blind=FALSE){
	if(requireNamespace("DESeq2")){
		# DESeq2 dataset
		o <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
																					colData = colData,
																					design = design)

		# Stabilize
		o <- DESeq2::varianceStabilizingTransformation(o, blind=blind)

		# Reshape
		o <- t(SummarizedExperiment::assay(o))

		# Return
		return(o)
	}else{
		stop("This function requires the DESeq2 package")
	}
}

#
# filter_and_normalize <- function(formula, data, counts, log=TRUE,
# 																 min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7, prior.count = 2, ...){
# 	# Setup up
# 	mod <- model.matrix(formula, data=data)
# 	o <- DGEList(counts)
#
# 	# Filter
# 	above <- filterByExpr(o,
# 												design=mod,
# 												min.count = min.count,
# 												min.total.count = min.total.count,
# 												large.n = large.n,
# 												min.prop = min.prop)
# 	o <- o[above,,keep.lib.sizes=FALSE]
#
# 	# Normalize
#
# 	# Log transform
# 	if(isTRUE(log)){
# 		o <- calcNormFactors(o, ...)
# 		o <- cpm(o, log=TRUE, prior.count = prior.count)
# 	}else{
# 		o <- o$counts
# 	}
#
# 	# tranpose
# 	o <- t(o)
#
# 	# Return
# 	o
# }
