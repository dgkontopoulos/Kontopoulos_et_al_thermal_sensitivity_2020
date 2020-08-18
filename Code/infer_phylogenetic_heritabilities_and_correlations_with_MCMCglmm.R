#!/usr/bin/Rscript

# This script estimates the phylogenetic heritability of each TPC parameter, 
# using MCMCglmm. It takes into account measurement uncertainty and the 
# covariances between pairs of TPC parameters. Missing values are 
# estimated based on the Missing-At-Random approach.
#
# It also prints the phenotypic, phylogenetically heritable, and residual 
# correlations between ln(E) and ln(W_op) for phytoplankton and 
# prokaryotes.

library(MCMCglmm)

#####################
# F U N C T I O N S #
#####################

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

# This function prepares TPC datasets (e.g., removes spaces in species 
# names, removes estimates from bad fits of the Sharpe-Schoolfield 
# model etc.) for the estimation of phylogenetic heritabilities. 
prepare_dataset <- function(dat, tree)
{
	dat$Species <- gsub("\\W", "_", dat$Species)
	
	dat <- dat[dat$R_squared >= 0.5 & dat$Species %in% tree$tip.label,]
	
	dat$B_0_fourth_root[dat$Data_points_rise < 4 | dat$B_0_fourth_root_error_variance == "NaN"] <- NA
	dat$B_0_fourth_root_error_variance[dat$Data_points_rise < 4 | dat$B_0_fourth_root_error_variance == "NaN"] <- NA
	
	dat$ln_E[dat$Data_points_rise < 4 | exp(dat$ln_E) > 4 | dat$ln_E_error_variance == "NaN"] <- NA
	dat$ln_E_error_variance[dat$Data_points_rise < 4 | exp(dat$ln_E) > 4 | dat$ln_E_error_variance == "NaN"] <- NA
	
	dat$ln_W_op[dat$Data_points_rise < 4 | dat$Data_points_fall < 2 | dat$ln_W_op_error_variance == "NaN"] <- NA
	dat$ln_W_op_error_variance[dat$Data_points_rise < 4 | dat$Data_points_fall < 2 | dat$ln_W_op_error_variance == "NaN"] <- NA
	
	dat$T_pk_squared[(dat$Data_points_rise < 2 | dat$Data_points_fall < 2) | dat$T_pk_squared_error_variance == "NaN"] <- NA
	dat$T_pk_squared_error_variance[(dat$Data_points_rise < 2 | dat$Data_points_fall < 2) | dat$T_pk_squared_error_variance == "NaN"] <- NA
	
	dat$ln_B_pk[dat$Data_points_rise < 2 | dat$Data_points_fall < 2 | dat$ln_B_pk_error_variance == "NaN"] <- NA
	dat$ln_B_pk_error_variance[dat$Data_points_rise < 2 | dat$Data_points_fall < 2 | dat$ln_B_pk_error_variance == "NaN"] <- NA
	
	dat$ln_E_D[dat$Data_points_fall < 4 | dat$ln_E_D_error_variance == "NaN"] <- NA
	dat$ln_E_D_error_variance[dat$Data_points_fall < 4 | dat$ln_E_D_error_variance == "NaN"] <- NA
	
	new_dat <- dat[, c(
			'Species', 'B_0_fourth_root', 'B_0_fourth_root_error_variance', 'ln_E', 
			'ln_E_error_variance', 'ln_W_op', 'ln_W_op_error_variance', 'T_pk_squared',
			'T_pk_squared_error_variance', 'ln_B_pk', 'ln_B_pk_error_variance', 'ln_E_D', 'ln_E_D_error_variance',
			'Data_points_rise', 'Data_points_fall'
		)
	]
	
	rows_to_drop <- c()
	for ( i in 1:nrow(new_dat) )
	{
		if ( 
			is.na(new_dat$B_0_fourth_root[i]) && 
			is.na(new_dat$ln_E[i]) &&
			is.na(new_dat$ln_W_op[i]) &&
			is.na(new_dat$T_pk_squared[i]) &&
			is.na(new_dat$ln_B_pk[i]) &&
			is.na(new_dat$ln_E_D[i])
		)
		{
			rows_to_drop <- c(rows_to_drop, i)
		}
	}
	
	if ( length(rows_to_drop) > 0 )
	{
		new_dat <- new_dat[-c(rows_to_drop),]
	}
	
	tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% unique(new_dat$Species))])
	
	return(list(dat = new_dat, tree = tree))
}

####################
# M A I N  C O D E #
####################

# Read the time-calibrated phylogeny.
tree <- read.nexus('../Data/final_calibrated_phylogeny.nex')
tree$node.label <- NULL

# Prepare a list of Cyanobacteria species in the datasets.
Cyanobacteria <- c(
	'Mastigocladus laminosus', 'Anabaena ucrainica',
	'Aphanizomenon flosaquae', 'Aphanizomenon gracile',
	'Aphanizomenon ovalisporum', 'Cylindrospermopsis raciborskii',
	'Limnothrix redekei', 'Microcystis aeruginosa',
	'Planktothrix agardhii', 'Prochlorococcus marinus',
	'Sphaerospermopsis aphanizomenoides', 'Spirulina platensis',
	'Synechococcus elongatus', 'Synechococcus lividus',
	'Trichodesmium erythraeum', 'Tychonema bourrellyi'
)

# Read the TPC datasets of phytoplankton and prokaryotes.
dat_phytoplankton <- read.csv('../Data/TPC_parameter_estimates_phytoplankton_r_max.csv', stringsAsFactors = FALSE)
dat_prokaryotes <- read.csv('../Data/TPC_parameter_estimates_prokaryotes_r_max.csv', stringsAsFactors = FALSE)

###################################################
# E U K AR Y O T I C    P H Y T O P L A N K T O N #
###################################################

# Remove Cyanobacteria from the phytoplankton dataset and estimate 
# phylogenetic heritabilities.
dat_eukaryotic_phytoplankton <- dat_phytoplankton[!(dat_phytoplankton$Species %in% Cyanobacteria),]
dat_eukaryotic_phytoplankton_cleaned <- prepare_dataset(dat_eukaryotic_phytoplankton, tree)

# Gather all the error variances together.
MEVs_eukaryotic_phytoplankton <- c(
	dat_eukaryotic_phytoplankton_cleaned$dat$B_0_fourth_root_error_variance, 
	dat_eukaryotic_phytoplankton_cleaned$dat$ln_E_error_variance, 
	dat_eukaryotic_phytoplankton_cleaned$dat$T_pk_squared_error_variance, 
	dat_eukaryotic_phytoplankton_cleaned$dat$ln_B_pk_error_variance, 
	dat_eukaryotic_phytoplankton_cleaned$dat$ln_E_D_error_variance, 
	dat_eukaryotic_phytoplankton_cleaned$dat$ln_W_op_error_variance
)
MEVs_eukaryotic_phytoplankton[is.na(MEVs_eukaryotic_phytoplankton)] <- 0

# Fit two chains of the intercept-only model.
set.seed(1)
model.1a_euk_phyto <- MCMCglmm(
	cbind(B_0_fourth_root, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait - 1,
	random=~us(trait):Species,
	family=rep("gaussian", 6),
	ginverse=list(Species=inverseA(dat_eukaryotic_phytoplankton_cleaned$tree, nodes = 'ALL', scale = TRUE)$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dat_eukaryotic_phytoplankton_cleaned$dat,
	mev = MEVs_eukaryotic_phytoplankton,
	rcov=~us(trait):units,
	nitt=200000000,
	burnin=20000000,
	thin=1000,
	verbose = FALSE
)

set.seed(2)
model.1b_euk_phyto <- MCMCglmm(
	cbind(B_0_fourth_root, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait - 1,
	random=~us(trait):Species,
	family=rep("gaussian", 6),
	ginverse=list(Species=inverseA(dat_eukaryotic_phytoplankton_cleaned$tree, nodes = 'ALL', scale = TRUE)$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dat_eukaryotic_phytoplankton_cleaned$dat,
	mev = MEVs_eukaryotic_phytoplankton,
	rcov=~us(trait):units,
	nitt=200000000,
	burnin=20000000,
	thin=1000,
	verbose = FALSE
)

# Ensure that the two chains have converged and that each parameter had 
# a sufficiently high effective sample size.
check_ESS(model.1a_euk_phyto, model.1b_euk_phyto)
check_convergence(model.1a_euk_phyto, model.1b_euk_phyto)

# Merge the estimates of the two chains.
model <- list(
	Sol = rbind(model.1a_euk_phyto$Sol, model.1b_euk_phyto$Sol),
	VCV = rbind(model.1a_euk_phyto$VCV, model.1b_euk_phyto$VCV)
)

# Calculate and report phylogenetic heritabilities.
herit_B_0_euk_phyto <- model$VCV[,"traitB_0_fourth_root:traitB_0_fourth_root.Species"] / 
	(model$VCV[,"traitB_0_fourth_root:traitB_0_fourth_root.Species"] + model$VCV[,"traitB_0_fourth_root:traitB_0_fourth_root.units"])
cat(
	"\n\nB_0 of eukaryotic phytoplankton: ", mean(herit_B_0_euk_phyto), " (", 
	HPDinterval(mcmc(herit_B_0_euk_phyto))[1], ",", 
	HPDinterval(mcmc(herit_B_0_euk_phyto))[2], ")\n", sep = ''
)

herit_E_euk_phyto <- model$VCV[,"traitln_E:traitln_E.Species"] / 
	(model$VCV[,"traitln_E:traitln_E.Species"] + model$VCV[,"traitln_E:traitln_E.units"])
cat(
	"E of eukaryotic phytoplankton: ", mean(herit_E_euk_phyto), " (", 
	HPDinterval(mcmc(herit_E_euk_phyto))[1], ",", 
	HPDinterval(mcmc(herit_E_euk_phyto))[2], ")\n", sep = ''
)

herit_T_pk_euk_phyto <- model$VCV[,"traitT_pk_squared:traitT_pk_squared.Species"] / 
	(model$VCV[,"traitT_pk_squared:traitT_pk_squared.Species"] + model$VCV[,"traitT_pk_squared:traitT_pk_squared.units"])
cat(
	"T_pk of eukaryotic phytoplankton: ", mean(herit_T_pk_euk_phyto), " (", 
	HPDinterval(mcmc(herit_T_pk_euk_phyto))[1], ",", 
	HPDinterval(mcmc(herit_T_pk_euk_phyto))[2], ")\n", sep = ''
)

herit_B_pk_euk_phyto <- model$VCV[,"traitln_B_pk:traitln_B_pk.Species"] / 
	(model$VCV[,"traitln_B_pk:traitln_B_pk.Species"] + model$VCV[,"traitln_B_pk:traitln_B_pk.units"])
cat(
	"B_pk of eukaryotic phytoplankton: ", mean(herit_B_pk_euk_phyto), " (", 
	HPDinterval(mcmc(herit_B_pk_euk_phyto))[1], ",", 
	HPDinterval(mcmc(herit_B_pk_euk_phyto))[2], ")\n", sep = ''
)

herit_E_D_euk_phyto <- model$VCV[,"traitln_E_D:traitln_E_D.Species"] / 
	(model$VCV[,"traitln_E_D:traitln_E_D.Species"] + model$VCV[,"traitln_E_D:traitln_E_D.units"])
cat(
	"E_D of eukaryotic phytoplankton: ", mean(herit_E_D_euk_phyto), " (", 
	HPDinterval(mcmc(herit_E_D_euk_phyto))[1], ",", 
	HPDinterval(mcmc(herit_E_D_euk_phyto))[2], ")\n", sep = ''
)

herit_W_op_euk_phyto <- model$VCV[,"traitln_W_op:traitln_W_op.Species"] / 
	(model$VCV[,"traitln_W_op:traitln_W_op.Species"] + model$VCV[,"traitln_W_op:traitln_W_op.units"])
cat(
	"W_op of eukaryotic phytoplankton: ", mean(herit_W_op_euk_phyto), " (", 
	HPDinterval(mcmc(herit_W_op_euk_phyto))[1], ",", 
	HPDinterval(mcmc(herit_W_op_euk_phyto))[2], ")\n\n\n", sep = ''
)

#############################
# P H Y T O P L A N K T O N #
#############################

# Put all Cyanobacteria in the phytoplankton dataset and estimate 
# phylogenetic heritabilities.
dat_phytoplankton <- rbind(dat_phytoplankton, dat_prokaryotes[dat_prokaryotes$Species %in% Cyanobacteria,])
dat_phytoplankton_cleaned <- prepare_dataset(dat_phytoplankton, tree)

# Gather all the error variances together.
MEVs_phytoplankton <- c(
	dat_phytoplankton_cleaned$dat$B_0_fourth_root_error_variance, 
	dat_phytoplankton_cleaned$dat$ln_E_error_variance, 
	dat_phytoplankton_cleaned$dat$T_pk_squared_error_variance, 
	dat_phytoplankton_cleaned$dat$ln_B_pk_error_variance, 
	dat_phytoplankton_cleaned$dat$ln_E_D_error_variance, 
	dat_phytoplankton_cleaned$dat$ln_W_op_error_variance
)
MEVs_phytoplankton[is.na(MEVs_phytoplankton)] <- 0

# Fit two chains of the intercept-only model.
set.seed(1)
model.1a_phyto <- MCMCglmm(
	cbind(B_0_fourth_root, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait - 1,
	random=~us(trait):Species,
	family=rep("gaussian", 6),
	ginverse=list(Species=inverseA(dat_phytoplankton_cleaned$tree, nodes = 'ALL', scale = TRUE)$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dat_phytoplankton_cleaned$dat,
	mev = MEVs_phytoplankton,
	rcov=~us(trait):units,
	nitt=200000000,
	burnin=20000000,
	thin=1000,
	verbose = FALSE
)

set.seed(2)
model.1b_phyto <- MCMCglmm(
	cbind(B_0_fourth_root, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait - 1,
	random=~us(trait):Species,
	family=rep("gaussian", 6),
	ginverse=list(Species=inverseA(dat_phytoplankton_cleaned$tree, nodes = 'ALL', scale = TRUE)$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dat_phytoplankton_cleaned$dat,
	mev = MEVs_phytoplankton,
	rcov=~us(trait):units,
	nitt=200000000,
	burnin=20000000,
	thin=1000,
	verbose = FALSE
)

# Ensure that the two chains have converged and that each parameter had 
# a sufficiently high effective sample size.
check_ESS(model.1a_phyto, model.1b_phyto)
check_convergence(model.1a_phyto, model.1b_phyto)

# Merge the estimates of the two chains.
model <- list(
	Sol = rbind(model.1a_phyto$Sol, model.1b_phyto$Sol),
	VCV = rbind(model.1a_phyto$VCV, model.1b_phyto$VCV)
)

# Calculate and report phylogenetic heritabilities.
herit_B_0_phyto <- model$VCV[,"traitB_0_fourth_root:traitB_0_fourth_root.Species"] / 
	(model$VCV[,"traitB_0_fourth_root:traitB_0_fourth_root.Species"] + model$VCV[,"traitB_0_fourth_root:traitB_0_fourth_root.units"])
cat(
	"\n\nB_0 of phytoplankton: ", mean(herit_B_0_phyto), " (", 
	HPDinterval(mcmc(herit_B_0_phyto))[1], ",", 
	HPDinterval(mcmc(herit_B_0_phyto))[2], ")\n", sep = ''
)

herit_E_phyto <- model$VCV[,"traitln_E:traitln_E.Species"] / 
	(model$VCV[,"traitln_E:traitln_E.Species"] + model$VCV[,"traitln_E:traitln_E.units"])
cat(
	"E of phytoplankton: ", mean(herit_E_phyto), " (", 
	HPDinterval(mcmc(herit_E_phyto))[1], ",", 
	HPDinterval(mcmc(herit_E_phyto))[2], ")\n", sep = ''
)

herit_T_pk_phyto <- model$VCV[,"traitT_pk_squared:traitT_pk_squared.Species"] / 
	(model$VCV[,"traitT_pk_squared:traitT_pk_squared.Species"] + model$VCV[,"traitT_pk_squared:traitT_pk_squared.units"])
cat(
	"T_pk of phytoplankton: ", mean(herit_T_pk_phyto), " (", 
	HPDinterval(mcmc(herit_T_pk_phyto))[1], ",", 
	HPDinterval(mcmc(herit_T_pk_phyto))[2], ")\n", sep = ''
)

herit_B_pk_phyto <- model$VCV[,"traitln_B_pk:traitln_B_pk.Species"] / 
	(model$VCV[,"traitln_B_pk:traitln_B_pk.Species"] + model$VCV[,"traitln_B_pk:traitln_B_pk.units"])
cat(
	"B_pk of phytoplankton: ", mean(herit_B_pk_phyto), " (", 
	HPDinterval(mcmc(herit_B_pk_phyto))[1], ",", 
	HPDinterval(mcmc(herit_B_pk_phyto))[2], ")\n", sep = ''
)

herit_E_D_phyto <- model$VCV[,"traitln_E_D:traitln_E_D.Species"] / 
	(model$VCV[,"traitln_E_D:traitln_E_D.Species"] + model$VCV[,"traitln_E_D:traitln_E_D.units"])
cat(
	"E_D of phytoplankton: ", mean(herit_E_D_phyto), " (", 
	HPDinterval(mcmc(herit_E_D_phyto))[1], ",", 
	HPDinterval(mcmc(herit_E_D_phyto))[2], ")\n", sep = ''
)

herit_W_op_phyto <- model$VCV[,"traitln_W_op:traitln_W_op.Species"] / 
	(model$VCV[,"traitln_W_op:traitln_W_op.Species"] + model$VCV[,"traitln_W_op:traitln_W_op.units"])
cat(
	"W_op of phytoplankton: ", mean(herit_W_op_phyto), " (", 
	HPDinterval(mcmc(herit_W_op_phyto))[1], ",", 
	HPDinterval(mcmc(herit_W_op_phyto))[2], ")\n\n\n", sep = ''
)

# Calculate and report correlations between ln(E) and ln(W_op).
E_W_op_cor_pheno_phyto <- (
	model$VCV[,"traitln_E:traitln_W_op.Species"] + model$VCV[,"traitln_E:traitln_W_op.units"]
)/sqrt(
	(
		model$VCV[,"traitln_E:traitln_E.Species"] + model$VCV[,"traitln_E:traitln_E.units"]
	) *
	(
		model$VCV[,"traitln_W_op:traitln_W_op.Species"] + model$VCV[,"traitln_W_op:traitln_W_op.units"]
	)
)
cat(
	"\nPhenotypic correlation between ln(E) and ln(W_op) among phytoplankton: ", mean(E_W_op_cor_pheno_phyto), " (", 
	HPDinterval(mcmc(E_W_op_cor_pheno_phyto))[1], ",", 
	HPDinterval(mcmc(E_W_op_cor_pheno_phyto))[2], ")\n", sep = ''
)

E_W_op_cor_phy_phyto <- (
	model$VCV[,"traitln_E:traitln_W_op.Species"]
)/sqrt(
	(
		model$VCV[,"traitln_E:traitln_E.Species"]
	) *
	(
		model$VCV[,"traitln_W_op:traitln_W_op.Species"]
	)
)
cat(
	"Phylogenetically heritable correlation between ln(E) and ln(W_op) among phytoplankton: ", mean(E_W_op_cor_phy_phyto), " (", 
	HPDinterval(mcmc(E_W_op_cor_phy_phyto))[1], ",", 
	HPDinterval(mcmc(E_W_op_cor_phy_phyto))[2], ")\n", sep = ''
)

E_W_op_cor_resid_phyto <- (
	model$VCV[,"traitln_E:traitln_W_op.units"]
)/sqrt(
	(
		model$VCV[,"traitln_E:traitln_E.units"]
	) *
	(
		model$VCV[,"traitln_W_op:traitln_W_op.units"]
	)
)
cat(
	"Residual correlation between ln(E) and ln(W_op) among phytoplankton: ", mean(E_W_op_cor_resid_phyto), " (", 
	HPDinterval(mcmc(E_W_op_cor_resid_phyto))[1], ",", 
	HPDinterval(mcmc(E_W_op_cor_resid_phyto))[2], ")\n\n\n", sep = ''
)

#########################
# P R O K A R Y O T E S #
#########################

# Remove all Cyanobacteria from the prokaryotes dataset and estimate 
# phylogenetic heritabilities.
dat_prokaryotes <- dat_prokaryotes[!(dat_prokaryotes$Species %in% Cyanobacteria),]
dat_prokaryotes_cleaned <- prepare_dataset(dat_prokaryotes, tree)

# Gather all the error variances together.
MEVs_prokaryotes <- c(
	dat_prokaryotes_cleaned$dat$B_0_fourth_root_error_variance, 
	dat_prokaryotes_cleaned$dat$ln_E_error_variance, 
	dat_prokaryotes_cleaned$dat$T_pk_squared_error_variance, 
	dat_prokaryotes_cleaned$dat$ln_B_pk_error_variance, 
	dat_prokaryotes_cleaned$dat$ln_E_D_error_variance, 
	dat_prokaryotes_cleaned$dat$ln_W_op_error_variance
)
MEVs_prokaryotes[is.na(MEVs_prokaryotes)] <- 0

# Fit two chains of the intercept-only model.
set.seed(1)
model.1a_prok <- MCMCglmm(
	cbind(B_0_fourth_root, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait - 1,
	random=~us(trait):Species,
	family=rep("gaussian", 6),
	ginverse=list(Species=inverseA(dat_prokaryotes_cleaned$tree, nodes = 'ALL', scale = TRUE)$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dat_prokaryotes_cleaned$dat,
	mev = MEVs_prokaryotes,
	rcov=~us(trait):units,
	nitt=200000000,
	burnin=20000000,
	thin=1000,
	verbose = FALSE
)

set.seed(2)
model.1b_prok <- MCMCglmm(
	cbind(B_0_fourth_root, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait - 1,
	random=~us(trait):Species,
	family=rep("gaussian", 6),
	ginverse=list(Species=inverseA(dat_prokaryotes_cleaned$tree, nodes = 'ALL', scale = TRUE)$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dat_prokaryotes_cleaned$dat,
	mev = MEVs_prokaryotes,
	rcov=~us(trait):units,
	nitt=200000000,
	burnin=20000000,
	thin=1000,
	verbose = FALSE
)

# Ensure that the two chains have converged and that each parameter had 
# a sufficiently high effective sample size.
check_ESS(model.1a_prok, model.1b_prok)
check_convergence(model.1a_prok, model.1b_prok)

# Merge the estimates of the two chains.
model <- list(
	Sol = rbind(model.1a_prok$Sol, model.1b_prok$Sol),
	VCV = rbind(model.1a_prok$VCV, model.1b_prok$VCV)
)

# Calculate and report phylogenetic heritabilities.
herit_B_0_prok <- model$VCV[,"traitB_0_fourth_root:traitB_0_fourth_root.Species"] / 
	(model$VCV[,"traitB_0_fourth_root:traitB_0_fourth_root.Species"] + model$VCV[,"traitB_0_fourth_root:traitB_0_fourth_root.units"])
cat(
	"\n\nB_0 of prokaryotes: ", mean(herit_B_0_prok), " (", 
	HPDinterval(mcmc(herit_B_0_prok))[1], ",", 
	HPDinterval(mcmc(herit_B_0_prok))[2], ")\n", sep = ''
)

herit_E_prok <- model$VCV[,"traitln_E:traitln_E.Species"] / 
	(model$VCV[,"traitln_E:traitln_E.Species"] + model$VCV[,"traitln_E:traitln_E.units"])
cat(
	"E of prokaryotes: ", mean(herit_E_prok), " (", 
	HPDinterval(mcmc(herit_E_prok))[1], ",", 
	HPDinterval(mcmc(herit_E_prok))[2], ")\n", sep = ''
)

herit_T_pk_prok <- model$VCV[,"traitT_pk_squared:traitT_pk_squared.Species"] / 
	(model$VCV[,"traitT_pk_squared:traitT_pk_squared.Species"] + model$VCV[,"traitT_pk_squared:traitT_pk_squared.units"])
cat(
	"T_pk of prokaryotes: ", mean(herit_T_pk_prok), " (", 
	HPDinterval(mcmc(herit_T_pk_prok))[1], ",", 
	HPDinterval(mcmc(herit_T_pk_prok))[2], ")\n", sep = ''
)

herit_B_pk_prok <- model$VCV[,"traitln_B_pk:traitln_B_pk.Species"] / 
	(model$VCV[,"traitln_B_pk:traitln_B_pk.Species"] + model$VCV[,"traitln_B_pk:traitln_B_pk.units"])
cat(
	"B_pk of prokaryotes: ", mean(herit_B_pk_prok), " (", 
	HPDinterval(mcmc(herit_B_pk_prok))[1], ",", 
	HPDinterval(mcmc(herit_B_pk_prok))[2], ")\n", sep = ''
)

herit_E_D_prok <- model$VCV[,"traitln_E_D:traitln_E_D.Species"] / 
	(model$VCV[,"traitln_E_D:traitln_E_D.Species"] + model$VCV[,"traitln_E_D:traitln_E_D.units"])
cat(
	"E_D of prokaryotes: ", mean(herit_E_D_prok), " (", 
	HPDinterval(mcmc(herit_E_D_prok))[1], ",", 
	HPDinterval(mcmc(herit_E_D_prok))[2], ")\n", sep = ''
)

herit_W_op_prok <- model$VCV[,"traitln_W_op:traitln_W_op.Species"] / 
	(model$VCV[,"traitln_W_op:traitln_W_op.Species"] + model$VCV[,"traitln_W_op:traitln_W_op.units"])
cat(
	"W_op of prokaryotes: ", mean(herit_W_op_prok), " (", 
	HPDinterval(mcmc(herit_W_op_prok))[1], ",", 
	HPDinterval(mcmc(herit_W_op_prok))[2], ")\n\n\n", sep = ''
)

# Calculate and report correlations between ln(E) and ln(W_op).
E_W_op_cor_pheno_prok <- (
	model$VCV[,"traitln_E:traitln_W_op.Species"] + model$VCV[,"traitln_E:traitln_W_op.units"]
)/sqrt(
	(
		model$VCV[,"traitln_E:traitln_E.Species"] + model$VCV[,"traitln_E:traitln_E.units"]
	) *
	(
		model$VCV[,"traitln_W_op:traitln_W_op.Species"] + model$VCV[,"traitln_W_op:traitln_W_op.units"]
	)
)
cat(
	"\nPhenotypic correlation between ln(E) and ln(W_op) among prokaryotes: ", mean(E_W_op_cor_pheno_prok), " (", 
	HPDinterval(mcmc(E_W_op_cor_pheno_prok))[1], ",", 
	HPDinterval(mcmc(E_W_op_cor_pheno_prok))[2], ")\n", sep = ''
)

E_W_op_cor_phy_prok <- (
	model$VCV[,"traitln_E:traitln_W_op.Species"]
)/sqrt(
	(
		model$VCV[,"traitln_E:traitln_E.Species"]
	) *
	(
		model$VCV[,"traitln_W_op:traitln_W_op.Species"]
	)
)
cat(
	"Phylogenetically heritable correlation between ln(E) and ln(W_op) among prokaryotes: ", mean(E_W_op_cor_phy_prok), " (", 
	HPDinterval(mcmc(E_W_op_cor_phy_prok))[1], ",", 
	HPDinterval(mcmc(E_W_op_cor_phy_prok))[2], ")\n", sep = ''
)

E_W_op_cor_resid_prok <- (
	model$VCV[,"traitln_E:traitln_W_op.units"]
)/sqrt(
	(
		model$VCV[,"traitln_E:traitln_E.units"]
	) *
	(
		model$VCV[,"traitln_W_op:traitln_W_op.units"]
	)
)
cat(
	"Residual correlation between ln(E) and ln(W_op) among prokaryotes: ", mean(E_W_op_cor_resid_prok), " (", 
	HPDinterval(mcmc(E_W_op_cor_resid_prok))[1], ",", 
	HPDinterval(mcmc(E_W_op_cor_resid_prok))[2], ")\n\n\n", sep = ''
)
