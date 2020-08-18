#!/usr/bin/Rscript

# This script fits phylogenetic and non-phylogenetic regressions of 
# ln(E) against binned latitude, trait identity, or just an intercept, 
# using MCMCglmm.
#
# 2 MCMC chains are run per model. The script checks that the chains have 
# converged on equivalent posterior distributions and that the 
# parameters are sufficiently sampled.
#
# The quality of the fitted models is evaluated according to the 
# Deviance Information Criterion (DIC). Marginal and conditional 
# R-squared values are also calculated.

library(MCMCglmm)

#####################
# F U N C T I O N S #
#####################

# This function returns the variance explained by the fixed effects in 
# an MCMCglmm.
calc_varF <- function(model)
{
	return(var(as.vector(apply(model$Sol,2,mean) %*% t(model$X))))
}

# This function calculates the marginal and conditional R-squared values 
# following Nakagawa & Schielzeth, Methods Ecol. Evol., 2013.
calc_r_squared <- function(model, varF)
{
	
	# Get the variance explained by random effects.
	varRandom <- mean(model$VCV[,'Species'])
	
	# Get the residual (unexplained) variance.
	varResid <- mean(model$VCV[,'units'])
	
	# Calculate the marginal and conditional R-squared values. 
	#
	# Marginal R-squared: fraction of variance explained by fixed effects only.
	# Conditional R-squared: fraction of variance explained by fixed and random effects.
	r_sq_marginal <- varF/(varF + varRandom + varResid)
	r_sq_conditional <- (varF + varRandom)/(varF + varRandom + varResid)
	
	return(list(marginal = r_sq_marginal, conditional = r_sq_conditional))
}

# This function checks if the two MCMCglmm chains converged on statistically 
# indistinguishable posterior distributions, based on the Potential 
# Scale Reduction Factor diagnostic (Gelman & Rubin, Stat. Sci., 1992).
check_convergence <- function(fit_1a, fit_1b)
{
	
	# Estimate the PSRF for each element of the fixed effects.
	# 
	# Values >= 1.1 indicate lack of convergence.
	# In that case, raise an error.
	gelman_Sol <- gelman.diag(mcmc.list(fit_1a$Sol, fit_1b$Sol), multivariate = FALSE)
	if ( length(gelman_Sol$psrf[
				!is.na(gelman_Sol$psrf[,1]) &
				gelman_Sol$psrf[,1] >= 1.1,
			]
		) > 0 
	)
	{
		stop("PROBLEM! The two runs have not converged!")
	}
	
	# Same thing for the random effects and the residual variance.
	gelman_VCV <- gelman.diag(mcmc.list(fit_1a$VCV, fit_1b$VCV), multivariate = FALSE)		
	if ( length(gelman_VCV$psrf[
				!is.na(gelman_VCV$psrf[,1]) &
				gelman_VCV$psrf[,1] >= 1.1,
			]
		) > 0 
	)
	{
		stop("PROBLEM! The two runs have not converged!")
	}
}

# This function checks that each parameter was adequately sampled from 
# the posterior, by calculating the Effective Sample Size.
check_ESS <- function(fit_1a, fit_1b)
{
	
	# Calculate the Effective Sample Size for each element of the 
	# fixed effects.
	merged_d_Sol <- rbind(fit_1a$Sol, fit_1b$Sol)
	ess_vals_Sol <- effectiveSize(merged_d_Sol)
	ess_vals_Sol <- ess_vals_Sol[ess_vals_Sol >= 1]
	
	# If there is a parameter with an Effective Sample Size below 
	# 200, raise an error.
	if ( length(ess_vals_Sol[ess_vals_Sol < 200]) > 0 )
	{
		stop("PROBLEM! The combined ESS is not big enough!")
	}

	# Same thing for the random effects and the residual variance.
	merged_d_VCV <- rbind(fit_1a$VCV, fit_1b$VCV)
	ess_vals_VCV <- effectiveSize(merged_d_VCV)
	ess_vals_VCV <- ess_vals_VCV[ess_vals_VCV >= 1]
	
	if ( length(ess_vals_VCV[ess_vals_VCV < 200]) > 0 )
	{
		stop("PROBLEM! The combined ESS is not big enough!")
	}	
}

# This function prepares the datasets (e.g., removes measurements for 
# which latitude of isolation is unknown) for the analysis.
prepare_E_vals <- function(dat, tree, trait)
{
	dat$Species <- gsub("\\W", "_", dat$Species)
	dat$Species <- gsub("_+$", "", dat$Species)
	
	dat <- dat[
		dat$R_squared >= 0.5 & 
		dat$Data_points_rise >= 4 & 
		exp(dat$ln_E) <= 4 & 
		dat$ln_E_error_variance != "NaN" &
		!(is.na(dat$Latitude)),
	]
	
	dat <- dat[dat$Species %in% tree$tip.label,]
	dat$Trait <- trait
	
	return(dat[,c("ln_E", "ln_E_error_variance", "Species", "Latitude", "Trait")])	
}


####################
# M A I N  C O D E #
####################

# Read the time-calibrated phylogeny.
tree <- read.nexus('../Data/final_calibrated_phylogeny.nex')
tree$node.label <- NULL

# Read the TPC datasets and prepare them for the analysis.
dat_photosynthesis <- read.csv('../Data/TPC_parameter_estimates_plants_net_photosynthesis_rate.csv', stringsAsFactors = FALSE)
dat_E_photosynthesis_cleaned <- prepare_E_vals(dat_photosynthesis, tree, 'net_photosynthesis_rate')

dat_respiration <- read.csv('../Data/TPC_parameter_estimates_plants_respiration_rate.csv', stringsAsFactors = FALSE)
dat_E_respiration_cleaned <- prepare_E_vals(dat_respiration, tree, 'respiration_rate')

dat_phytoplankton <- read.csv('../Data/TPC_parameter_estimates_phytoplankton_r_max.csv', stringsAsFactors = FALSE)
dat_E_phytoplankton_cleaned <- prepare_E_vals(dat_phytoplankton, tree, 'r_max')

dat_prokaryotes <- read.csv('../Data/TPC_parameter_estimates_prokaryotes_r_max.csv', stringsAsFactors = FALSE)
dat_E_prokaryotes_cleaned <- prepare_E_vals(dat_prokaryotes, tree, 'r_max')

# Merge the TPC dataset.
all_dat_E <- rbind(
	dat_E_photosynthesis_cleaned,	dat_E_respiration_cleaned, 
	dat_E_prokaryotes_cleaned,		dat_E_phytoplankton_cleaned
)

# Set r_max as the baseline trait for comparisons.
all_dat_E$Trait <- factor(
	all_dat_E$Trait, levels = c(
		'r_max', 'net_photosynthesis_rate', 'respiration_rate'
	)
)

# Prepare the latitude variable.
all_dat_E$Latitude_bins <- all_dat_E$Latitude
all_dat_E$Latitude_bins[all_dat_E$Latitude >= 0 & all_dat_E$Latitude < 30] <- '0-30'
all_dat_E$Latitude_bins[all_dat_E$Latitude >= 30 & all_dat_E$Latitude < 60] <- '30-60'
all_dat_E$Latitude_bins[all_dat_E$Latitude >= 60] <- '>60'
all_dat_E$Latitude_bins[all_dat_E$Latitude < 0 & all_dat_E$Latitude > -30] <- '0-30'
all_dat_E$Latitude_bins[all_dat_E$Latitude <= -30 & all_dat_E$Latitude > -60] <- '30-60'
all_dat_E$Latitude_bins[all_dat_E$Latitude <= -60] <- '>60'

all_dat_E$Latitude_bins <- factor(
	all_dat_E$Latitude_bins, levels = c('30-60', '0-30', '>60')
)

# Drop any species from the phylogeny for which there are no E estimates.
tree_E <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% unique(all_dat_E$Species))])

# Fit the models.
set.seed(7)
fit_4a <- MCMCglmm(ln_E ~ Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree_E, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)

set.seed(8)
fit_4b <- MCMCglmm(ln_E ~ Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree_E, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)
check_ESS(fit_4a, fit_4b)
check_convergence(fit_4a, fit_4b)

set.seed(9)
fit_5a <- MCMCglmm(ln_E ~ 1, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree_E, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)

set.seed(10)
fit_5b <- MCMCglmm(ln_E ~ 1, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree_E, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)
check_ESS(fit_5a, fit_5b)
check_convergence(fit_5a, fit_5b)

set.seed(11)
fit_6a <- MCMCglmm(ln_E ~ Latitude_bins + Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree_E, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)

set.seed(12)
fit_6b <- MCMCglmm(ln_E ~ Latitude_bins + Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree_E, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)
check_ESS(fit_6a, fit_6b)
check_convergence(fit_6a, fit_6b)

set.seed(13)
fit_7a <- MCMCglmm(ln_E ~ Latitude_bins*Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree_E, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)

set.seed(14)
fit_7b <- MCMCglmm(ln_E ~ Latitude_bins*Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree_E, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)
check_ESS(fit_7a, fit_7b)
check_convergence(fit_7a, fit_7b)

set.seed(15)
fit_8a <- MCMCglmm(ln_E ~ Latitude_bins, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree_E, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)

set.seed(16)
fit_8b <- MCMCglmm(ln_E ~ Latitude_bins, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree_E, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)
check_ESS(fit_8a, fit_8b)
check_convergence(fit_8a, fit_8b)

set.seed(23)
fit_4a_NON_PHYLO <- MCMCglmm(ln_E ~ Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)

set.seed(24)
fit_4b_NON_PHYLO <- MCMCglmm(ln_E ~ Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)
check_ESS(fit_4a_NON_PHYLO, fit_4b_NON_PHYLO)
check_convergence(fit_4a_NON_PHYLO, fit_4b_NON_PHYLO)

set.seed(25)
fit_5a_NON_PHYLO <- MCMCglmm(ln_E ~ 1, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)

set.seed(26)
fit_5b_NON_PHYLO <- MCMCglmm(ln_E ~ 1, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)
check_ESS(fit_5a_NON_PHYLO, fit_5b_NON_PHYLO)
check_convergence(fit_5a_NON_PHYLO, fit_5b_NON_PHYLO)

set.seed(27)
fit_6a_NON_PHYLO <- MCMCglmm(ln_E ~ Latitude_bins + Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)

set.seed(28)
fit_6b_NON_PHYLO <- MCMCglmm(ln_E ~ Latitude_bins + Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)
check_ESS(fit_6a_NON_PHYLO, fit_6b_NON_PHYLO)
check_convergence(fit_6a_NON_PHYLO, fit_6b_NON_PHYLO)

set.seed(29)
fit_7a_NON_PHYLO <- MCMCglmm(ln_E ~ Latitude_bins*Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)

set.seed(30)
fit_7b_NON_PHYLO <- MCMCglmm(ln_E ~ Latitude_bins*Trait, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)
check_ESS(fit_7a_NON_PHYLO, fit_7b_NON_PHYLO)
check_convergence(fit_7a_NON_PHYLO, fit_7b_NON_PHYLO)

set.seed(31)
fit_8a_NON_PHYLO <- MCMCglmm(ln_E ~ Latitude_bins, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)

set.seed(32)
fit_8b_NON_PHYLO <- MCMCglmm(ln_E ~ Latitude_bins, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = all_dat_E,
    nitt = 5000000,
    burnin = 500000,
    thin = 1000,
    verbose = TRUE,
    mev = all_dat_E$ln_E_error_variance
)
check_ESS(fit_8a_NON_PHYLO, fit_8b_NON_PHYLO)
check_convergence(fit_8a_NON_PHYLO, fit_8b_NON_PHYLO)

# The best-fitting model is the 6th non-phylogenetic one, i.e., with 
# fixed effects of binned latitude and trait identiy, and with a 
# species random effect on the intercept.

# Merge the estimates of the two chains belonging to the best model.
joined_best_fit <- list(
	Sol = rbind(fit_6a_NON_PHYLO$Sol, fit_6b_NON_PHYLO$Sol),
	VCV = rbind(fit_6a_NON_PHYLO$VCV, fit_6b_NON_PHYLO$VCV)
)

# Calculate the marginal and conditional R-squared values.
r_squared_vals <- calc_r_squared(
	joined_best_fit, mean(c(calc_varF(fit_6a_NON_PHYLO), calc_varF(fit_6b_NON_PHYLO)))
)
