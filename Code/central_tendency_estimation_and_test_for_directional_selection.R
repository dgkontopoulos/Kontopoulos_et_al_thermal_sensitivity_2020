#!/usr/bin/Rscript

# This script fits Eq. 2 from the main text, i.e., it estimates the 
# central tendency of thermal sensitivity and also quantifies a 
# potential trend towards lower or higher values.
#
# This script requires the output files produced from fitting the stable 
# model. Thus, to successfully run this script, make sure that you run 
# the fit_stable_model.R script first.

library(geiger)
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

# This function fits models to quantify the central tendency of E 
# and a putative trend towards lower or higher values with time.
fit_models_E <- function(tree, dat, working_dir)
{
	
	# Prepare the dataset (e.g., remove spaces in species names, 
	# remove estimates from bad fits of the Sharpe-Schoolfield model etc.).
	dat$Species <- gsub("\\W", "_", dat$Species)
	dat$Species <- gsub("_+$", "", dat$Species)
	
	dat <- dat[dat$R_squared >= 0.5 & dat$Data_points_rise >= 4 & exp(dat$ln_E) <= 4 & dat$ln_E_error_variance != "NaN",]
	
	dat <- dat[dat$Species %in% tree$tip.label,]
	tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% unique(dat$Species))])

	ids_to_keep <- c()
	for ( i in unique(dat$Species) )
	{
		temp_dat <- dat[dat$Species == i,]
		
		# If there are multiple fits for the same species, keep the 
		# one with the highest R-squared value.
		new_id <- row.names(temp_dat[temp_dat$R_squared == max(temp_dat$R_squared),])
		if ( length(new_id) > 1 )
		{
			new_id <- new_id[1]
		}
		
		ids_to_keep <- c(
			ids_to_keep, 
			new_id
		)
	}
	
	dat <- dat[as.character(ids_to_keep),]
	
	# Read the stable model fits.
	chain1 <- read.delim(paste(working_dir, 'out.chain1.log', sep = ''))
	chain2 <- read.delim(paste(working_dir, 'out.chain2.log', sep = ''))
	chain3 <- read.delim(paste(working_dir, 'out.chain3.log', sep = ''))
	chain4 <- read.delim(paste(working_dir, 'out.chain4.log', sep = ''))
	
	chain1 <- chain1[(round(nrow(chain1) * 0.25) + 1):nrow(chain1),2:ncol(chain1)]
	chain2 <- chain2[(round(nrow(chain2) * 0.25) + 1):nrow(chain2),2:ncol(chain2)]
	chain3 <- chain3[(round(nrow(chain3) * 0.25) + 1):nrow(chain3),2:ncol(chain3)]
	chain4 <- chain4[(round(nrow(chain4) * 0.25) + 1):nrow(chain4),2:ncol(chain4)]
	
	all_samples <- rbind(chain1, chain2, chain3, chain4)
	ancestral_states <- colMeans(all_samples[,4:(ncol(all_samples) - 1)])
	names(ancestral_states) <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)
	
	sp_estimates <- dat$ln_E
	names(sp_estimates) <- dat$Species
	
	tree$node.label <- names(ancestral_states)

	# Collect node ages. 
	temp_dataset <- data.frame(
		all_estimates = c(sp_estimates, ancestral_states),
		time = c(rep(0, length(sp_estimates)), 1 - branching.times(rescale(tree, 'depth', 1))),
		Node = c(names(sp_estimates), names(ancestral_states))
	)
	temp_dataset <- temp_dataset[temp_dataset$Node != names(ancestral_states)[1],]
	
	# Fit phylogenetic and non-phylogenetic models.
	set.seed(1)
	fit_PHYLO_a <- MCMCglmm(all_estimates ~ time, family = "gaussian",
	    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	    random=~Node,
	    data = temp_dataset,
	    nitt = 1000000,
	    burnin = 100000,
	    thin = 1000,
	    ginverse=list(Node=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
	    verbose = TRUE
	)
	set.seed(2)
	fit_PHYLO_b <- MCMCglmm(all_estimates ~ time, family = "gaussian",
	    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	    random=~Node,
	    data = temp_dataset,
	    nitt = 1000000,
	    burnin = 100000,
	    thin = 1000,
	    ginverse=list(Node=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
	    verbose = TRUE
	)
	check_ESS(fit_PHYLO_a, fit_PHYLO_b)
	check_convergence(fit_PHYLO_a, fit_PHYLO_b)
	
	set.seed(3)
	fit_NON_PHYLO_a <- MCMCglmm(all_estimates ~ time, family = "gaussian",
	    prior=list(R=list(V=diag(1),nu=1.002)),
	    data = temp_dataset,
	    nitt = 1000000,
	    burnin = 100000,
	    thin = 1000,
	    verbose = TRUE
	)
	set.seed(4)
	fit_NON_PHYLO_b <- MCMCglmm(all_estimates ~ time, family = "gaussian",
	    prior=list(R=list(V=diag(1),nu=1.002)),
	    data = temp_dataset,
	    nitt = 1000000,
	    burnin = 100000,
	    thin = 1000,
	    verbose = TRUE
	)
	check_ESS(fit_NON_PHYLO_a, fit_NON_PHYLO_b)
	check_convergence(fit_NON_PHYLO_a, fit_NON_PHYLO_b)
	
	all_PHYLO <- list(
		Sol = rbind(fit_PHYLO_a$Sol, fit_PHYLO_b$Sol),
		VCV = rbind(fit_PHYLO_a$VCV, fit_PHYLO_b$VCV)
	)
	all_NON_PHYLO <- list(
		Sol = rbind(fit_NON_PHYLO_a$Sol, fit_NON_PHYLO_b$Sol),
		VCV = rbind(fit_NON_PHYLO_a$VCV, fit_NON_PHYLO_b$VCV)
	)
	
	# Save the model objects to a file, identify the best-fitting models
	# on the basis of DIC, and summarise the results.
	save(
		all_PHYLO, 
		all_NON_PHYLO, 
		file = paste(working_dir, 'objs.Rda', sep = '')
	)
	
	sink(file = paste(working_dir, 'summary.txt', sep = ''))
	
	DIC_phylo <- mean(c(fit_PHYLO_a$DIC, fit_PHYLO_b$DIC))
	DIC_non_phylo <- mean(c(fit_NON_PHYLO_a$DIC, fit_NON_PHYLO_b$DIC))
	
	cat("Mean DIC of phylogenetic model: ", DIC_phylo, "\n")
	cat("Mean DIC of non-phylogenetic model: ", DIC_non_phylo, "\n\n")
	
	if ( DIC_phylo < DIC_non_phylo )
	{
		print(colMeans(all_PHYLO$Sol))
		cat("\n")
		print(HPDinterval(mcmc(all_PHYLO$Sol)))
		cat("\n")
		print(colMeans(all_PHYLO$VCV))
		cat("\n")
		print(HPDinterval(mcmc(all_PHYLO$VCV)))
	} else
	{
		print(colMeans(all_NON_PHYLO$Sol))
		cat("\n")
		print(HPDinterval(mcmc(all_NON_PHYLO$Sol)))
		cat("\n")
		print(colMeans(all_NON_PHYLO$VCV))
		cat("\n")
		print(HPDinterval(mcmc(all_NON_PHYLO$VCV)))		
	}
	
	sink()
	return()
}

# This function fits models to quantify the central tendency of W_op 
# and a putative trend towards lower or higher values with time.
fit_models_W_op <- function(tree, dat, working_dir)
{
	
	# Prepare the dataset (e.g., remove spaces in species names, 
	# remove estimates from bad fits of the Sharpe-Schoolfield model etc.).
	dat$Species <- gsub("\\W", "_", dat$Species)
	dat$Species <- gsub("_+$", "", dat$Species)
	
	dat <- dat[dat$R_squared >= 0.5 & dat$Data_points_rise >= 4 & dat$Data_points_fall >= 2 & dat$ln_W_op_error_variance != "NaN",]
	
	dat <- dat[dat$Species %in% tree$tip.label,]
	tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% unique(dat$Species))])

	ids_to_keep <- c()
	for ( i in unique(dat$Species) )
	{
		temp_dat <- dat[dat$Species == i,]
		
		# If there are multiple fits for the same species, keep the 
		# one with the highest R-squared value.
		new_id <- row.names(temp_dat[temp_dat$R_squared == max(temp_dat$R_squared),])
		if ( length(new_id) > 1 )
		{
			new_id <- new_id[1]
		}
		
		ids_to_keep <- c(
			ids_to_keep, 
			new_id
		)
	}
	
	dat <- dat[as.character(ids_to_keep),]

	# Read the stable model fits.	
	chain1 <- read.delim(paste(working_dir, 'out.chain1.log', sep = ''))
	chain2 <- read.delim(paste(working_dir, 'out.chain2.log', sep = ''))
	chain3 <- read.delim(paste(working_dir, 'out.chain3.log', sep = ''))
	chain4 <- read.delim(paste(working_dir, 'out.chain4.log', sep = ''))
	
	chain1 <- chain1[(round(nrow(chain1) * 0.25) + 1):nrow(chain1),2:ncol(chain1)]
	chain2 <- chain2[(round(nrow(chain2) * 0.25) + 1):nrow(chain2),2:ncol(chain2)]
	chain3 <- chain3[(round(nrow(chain3) * 0.25) + 1):nrow(chain3),2:ncol(chain3)]
	chain4 <- chain4[(round(nrow(chain4) * 0.25) + 1):nrow(chain4),2:ncol(chain4)]
	
	all_samples <- rbind(chain1, chain2, chain3, chain4)
	ancestral_states <- colMeans(all_samples[,4:(ncol(all_samples) - 1)])
	names(ancestral_states) <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)
	
	sp_estimates <- dat$ln_W_op
	names(sp_estimates) <- dat$Species
	
	tree$node.label <- names(ancestral_states)
	
	# Collect node ages.
	temp_dataset <- data.frame(
		all_estimates = c(sp_estimates, ancestral_states),
		time = c(rep(0, length(sp_estimates)), 1 - branching.times(rescale(tree, 'depth', 1))),
		Node = c(names(sp_estimates), names(ancestral_states))
	)
	temp_dataset <- temp_dataset[temp_dataset$Node != names(ancestral_states)[1],]
	
	# Fit phylogenetic and non-phylogenetic models.
	set.seed(1)
	fit_PHYLO_a <- MCMCglmm(all_estimates ~ time, family = "gaussian",
	    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	    random=~Node,
	    data = temp_dataset,
	    nitt = 1000000,
	    burnin = 100000,
	    thin = 1000,
	    ginverse=list(Node=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
	    verbose = TRUE
	)
	set.seed(2)
	fit_PHYLO_b <- MCMCglmm(all_estimates ~ time, family = "gaussian",
	    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	    random=~Node,
	    data = temp_dataset,
	    nitt = 1000000,
	    burnin = 100000,
	    thin = 1000,
	    ginverse=list(Node=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
	    verbose = TRUE
	)
	check_ESS(fit_PHYLO_a, fit_PHYLO_b)
	check_convergence(fit_PHYLO_a, fit_PHYLO_b)
	
	set.seed(3)
	fit_NON_PHYLO_a <- MCMCglmm(all_estimates ~ time, family = "gaussian",
	    prior=list(R=list(V=diag(1),nu=1.002)),
	    data = temp_dataset,
	    nitt = 1000000,
	    burnin = 100000,
	    thin = 1000,
	    verbose = TRUE
	)
	set.seed(4)
	fit_NON_PHYLO_b <- MCMCglmm(all_estimates ~ time, family = "gaussian",
	    prior=list(R=list(V=diag(1),nu=1.002)),
	    data = temp_dataset,
	    nitt = 1000000,
	    burnin = 100000,
	    thin = 1000,
	    verbose = TRUE
	)
	check_ESS(fit_NON_PHYLO_a, fit_NON_PHYLO_b)
	check_convergence(fit_NON_PHYLO_a, fit_NON_PHYLO_b)
	
	all_PHYLO <- list(
		Sol = rbind(fit_PHYLO_a$Sol, fit_PHYLO_b$Sol),
		VCV = rbind(fit_PHYLO_a$VCV, fit_PHYLO_b$VCV)
	)
	all_NON_PHYLO <- list(
		Sol = rbind(fit_NON_PHYLO_a$Sol, fit_NON_PHYLO_b$Sol),
		VCV = rbind(fit_NON_PHYLO_a$VCV, fit_NON_PHYLO_b$VCV)
	)
	
	save(
		all_PHYLO, 
		all_NON_PHYLO, 
		file = paste(working_dir, 'objs.Rda', sep = '')
	)
	
	# Save the model objects to a file, identify the best-fitting models
	# on the basis of DIC, and summarise the results.
	sink(file = paste(working_dir, 'summary.txt', sep = ''))
	
	DIC_phylo <- mean(c(fit_PHYLO_a$DIC, fit_PHYLO_b$DIC))
	DIC_non_phylo <- mean(c(fit_NON_PHYLO_a$DIC, fit_NON_PHYLO_b$DIC))
	
	cat("Mean DIC of phylogenetic model: ", DIC_phylo, "\n")
	cat("Mean DIC of non-phylogenetic model: ", DIC_non_phylo, "\n\n")
	
	if ( DIC_phylo < DIC_non_phylo )
	{
		print(colMeans(all_PHYLO$Sol))
		cat("\n")
		print(HPDinterval(mcmc(all_PHYLO$Sol)))
		cat("\n")
		print(colMeans(all_PHYLO$VCV))
		cat("\n")
		print(HPDinterval(mcmc(all_PHYLO$VCV)))
	} else
	{
		print(colMeans(all_NON_PHYLO$Sol))
		cat("\n")
		print(HPDinterval(mcmc(all_NON_PHYLO$Sol)))
		cat("\n")
		print(colMeans(all_NON_PHYLO$VCV))
		cat("\n")
		print(HPDinterval(mcmc(all_NON_PHYLO$VCV)))		
	}
	
	sink()
	return()
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

# Read the TPC datasets of phytoplankton, prokaryotes, and plants.
dat_phytoplankton <- read.csv('../Data/TPC_parameter_estimates_phytoplankton_r_max.csv', stringsAsFactors = FALSE)
dat_prokaryotes <- read.csv('../Data/TPC_parameter_estimates_prokaryotes_r_max.csv', stringsAsFactors = FALSE)
dat_photosynthesis <- read.csv('../Data/TPC_parameter_estimates_plants_net_photosynthesis_rate.csv', stringsAsFactors = FALSE)
dat_respiration <- read.csv('../Data/TPC_parameter_estimates_plants_respiration_rate.csv', stringsAsFactors = FALSE)

# Put all Cyanobacteria in the phytoplankton dataset and fit the 
# models.
dat_phytoplankton <- rbind(dat_phytoplankton, dat_prokaryotes[dat_prokaryotes$Species %in% Cyanobacteria,])
fit_models_E(
	tree, dat_phytoplankton, '../Results/stable/E/phytoplankton/'
)
fit_models_W_op(
	tree, dat_phytoplankton, '../Results/stable/W_op/phytoplankton/'
)

# Remove all Cyanobacteria from the prokaryotes dataset and fit the 
# models.
dat_prokaryotes <- dat_prokaryotes[!(dat_prokaryotes$Species %in% Cyanobacteria),]
fit_models_E(
	tree, dat_prokaryotes, '../Results/stable/E/prokaryotes/'
)
fit_models_W_op(
	tree, dat_prokaryotes, '../Results/stable/W_op/prokaryotes/'
)

# Fit the models also to the photosynthesis and respiration datasets.
fit_models_E(
	tree, dat_photosynthesis, '../Results/stable/E/photosynthesis/'
)
fit_models_W_op(
	tree, dat_photosynthesis, '../Results/stable/W_op/photosynthesis/'
)

fit_models_E(
	tree, dat_respiration, '../Results/stable/E/respiration/'
)
fit_models_W_op(
	tree, dat_respiration, '../Results/stable/W_op/respiration/'
)
