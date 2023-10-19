#### Custom strain function ####
# calcStrain <- function(object){
# 	# Rank
# 	r <- object$CCA$rank
#
# 	# Extract scores
# 	o <- scores(object, display=c("wa", "lc"), scaling = "sites", choices = seq_len(r))
#
# 	# Calc strain
# 	o <- vapply(1:nrow(o$sites),
# 							function(i) dist(rbind(o$sites[i,], o$constraints[i,]))[1],
# 							numeric(1))
#
# 	# Return
# 	o
# }

#### Built-in sample quality ####

# calcSampleQuality <- function(object){
# 	data.frame(Constrained=inertcomp(object, display="sites", unity=FALSE, proportional = FALSE)[,"CCA"],
# 						 Unity=inertcomp(object, display="sites", unity = TRUE, proportional = FALSE)[,"CCA"],
# 						 Proportional=inertcomp(object, display="sites", unity=FALSE, proportional = TRUE)[,"CCA"],
# 						 GoF=goodness(object, display="sites", summarize = TRUE),
# 						 Strain=calcStrain(object))
# }
# calcSampleQuality <- function(object){
# 	data.frame(Inertia=inertcomp(object, display="sites", unity=FALSE, proportional = FALSE)[,"CCA"], # unity is just a scaling
# 						 Proportional=inertcomp(object, display="sites", unity=FALSE, proportional = TRUE)[,"CCA"]) # same as goodness
# 						 #Strain=calcStrain(object))
# }
