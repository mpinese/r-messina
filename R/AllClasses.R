setOldClass("Surv")

setClassUnion("MessinaY", c("Surv", "vector"))

.MessinaParameters <- setClass(	"MessinaParameters",
								slots = c(	x = "matrix",
											y = "MessinaY",
											features = "vector",
											samples = "vector",
											perf_requirement = "list",		# For MessinaSurvResult, list with names objective_type, min_objective.  For MessinaClassResult, list with names min_sensitivity, min_specificity
											minimum_group_fraction = "numeric",
											training_fraction = "numeric",
											num_bootstraps = "integer",
											prng_seed = "integer"))

.MessinaFits <- setClass("MessinaFits",
						slots = c(	summary = "data.frame",					# columns: passed, type, threshold, posk, margin, ptrue, psuccessful
									objective_surfaces = "list"))			# list of data.frames, each with columns: cutoff, objective (MessinaSurvResult), or cutoff, sensitivity, specificity (MessinaClassResult)

.MessinaResult <- setClass(	"MessinaResult", 
							slots = c(	problem_type = "character",
										parameters = "MessinaParameters",
										perf_estimates = "data.frame",		# Columns for MessinaClassResult: mean_tpr, mean_fpr, mean_tnr, mean_fnr, var_tpr, var_fpr, var_tnr, var_fnr, mean_sens, mean_spec
																			# Columns for MessinaSurvResult: mean_obj
										fits = "MessinaFits"))

.MessinaClassResult <- setClass("MessinaClassResult", contains = "MessinaResult")

.MessinaSurvResult <- setClass(	"MessinaSurvResult", contains = "MessinaResult")
