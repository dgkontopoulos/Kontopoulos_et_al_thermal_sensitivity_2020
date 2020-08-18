#!/usr/bin/Rscript

# This script estimates the phylogenetic heritability (Pagel's lambda) 
# separately for each TPC parameter of phytoplankton and prokaryotes, 
# using the BayesTraits program.

library(ape)
library(coda)

#####################
# F U N C T I O N S #
#####################

# This function checks if the two MCMCglmm chains converged on statistically 
# indistinguishable posterior distributions, based on the Potential 
# Scale Reduction Factor diagnostic (Gelman & Rubin, Stat. Sci., 1992).
check_convergence <- function(chain1, chain2)
{
	
	# Estimate the PSRF for each parameter.
	# 
	# Values >= 1.1 indicate lack of convergence.
	# In that case, raise an error.
	gelman_result <- gelman.diag(mcmc.list(mcmc(chain1), mcmc(chain2)), multivariate = FALSE)
	if ( length(gelman_result$psrf[
				!is.na(gelman_result$psrf[,1]) &
				gelman_result$psrf[,1] >= 1.1,
			]
		) > 0 
	)
	{
		stop("PROBLEM! The two chains have not converged!")
	}
}

# This function checks that each parameter was adequately sampled from 
# the posterior, by calculating the Effective Sample Size.
check_ESS <- function(chain1, chain2)
{

	# Calculate the Effective Sample Size for each parameter.
	ess_chain1 <- effectiveSize(chain1)
	ess_chain1 <- ess_chain1[ess_chain1 > 0]
	
	ess_chain2 <- effectiveSize(chain2)
	ess_chain2 <- ess_chain2[ess_chain2 > 0]

	# If there is a parameter with an Effective Sample Size below 
	# 200, raise an error.
	if ( min(ess_chain1) < 200 || min(ess_chain2) < 200 )
	{
		stop("PROBLEM! The ESS is not big enough!")
	}	
}

# This function prepares TPC datasets (e.g., removes spaces in species 
# names, removes estimates from bad fits of the Sharpe-Schoolfield 
# model etc.) for the estimation of Pagel's lambda with BayesTraits.
prepare_dataset <- function(dat, tree, working_dir)
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
	
	# Create a directory (if not already available) to store the 
	# input and output files.
	dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
	
	# Prepare the datasets required for BayesTraits to run the analysis.
	df_for_BayesTraits <- data.frame(
		Species = rep('-', length(unique(new_dat$Species))),
		Linkage = rep('Unlinked', length(unique(new_dat$Species))),
		B_0_fourth_root = rep('-', length(unique(new_dat$Species))),
		ln_E = rep('-', length(unique(new_dat$Species))),
		ln_W_op = rep('-', length(unique(new_dat$Species))),
		T_pk_squared = rep('-', length(unique(new_dat$Species))),
		ln_B_pk = rep('-', length(unique(new_dat$Species))),
		ln_E_D = rep('-', length(unique(new_dat$Species)))
	)
	
	df_for_BayesTraits$Species <- as.character(df_for_BayesTraits$Species)
	df_for_BayesTraits$Linkage <- as.character(df_for_BayesTraits$Linkage)
	df_for_BayesTraits$B_0_fourth_root <- as.character(df_for_BayesTraits$B_0_fourth_root)
	df_for_BayesTraits$ln_E <- as.character(df_for_BayesTraits$ln_E)
	df_for_BayesTraits$ln_W_op <- as.character(df_for_BayesTraits$ln_W_op)
	df_for_BayesTraits$T_pk_squared <- as.character(df_for_BayesTraits$T_pk_squared)
	df_for_BayesTraits$ln_B_pk <- as.character(df_for_BayesTraits$ln_B_pk)
	df_for_BayesTraits$ln_E_D <- as.character(df_for_BayesTraits$ln_E_D)

	df_for_BayesTraits2 <- data.frame(
		Species = rep('?', length(unique(new_dat$Species))),
		B_0_fourth_root = rep('?', length(unique(new_dat$Species))),
		ln_E = rep('?', length(unique(new_dat$Species))),
		ln_W_op = rep('?', length(unique(new_dat$Species))),
		T_pk_squared = rep('?', length(unique(new_dat$Species))),
		ln_B_pk = rep('?', length(unique(new_dat$Species))),
		ln_E_D = rep('?', length(unique(new_dat$Species)))
	)
	
	# Prepare a data frame to hold multiple estimates per TPC parameter 
	# per species, when available.
	df_for_BayesTraits2$Species <- as.character(df_for_BayesTraits2$Species)
	df_for_BayesTraits2$B_0_fourth_root <- as.character(df_for_BayesTraits2$B_0_fourth_root)
	df_for_BayesTraits2$ln_E <- as.character(df_for_BayesTraits2$ln_E)
	df_for_BayesTraits2$ln_W_op <- as.character(df_for_BayesTraits2$ln_W_op)
	df_for_BayesTraits2$T_pk_squared <- as.character(df_for_BayesTraits2$T_pk_squared)
	df_for_BayesTraits2$ln_B_pk <- as.character(df_for_BayesTraits2$ln_B_pk)
	df_for_BayesTraits2$ln_E_D <- as.character(df_for_BayesTraits2$ln_E_D)
	
	counter <- 1
	for ( i in unique(new_dat$Species) )
	{
		temp_new_dat <- new_dat[new_dat$Species == i,]
		
		df_for_BayesTraits$Species[counter] <- i
		df_for_BayesTraits2$Species[counter] <- i
		
		if ( length(na.omit(temp_new_dat$B_0_fourth_root)) > 0 )
		{
			df_for_BayesTraits$B_0_fourth_root[counter] <- paste(
				na.omit(temp_new_dat$B_0_fourth_root), collapse = ','
			)
			df_for_BayesTraits2$B_0_fourth_root[counter] <- as.character(na.omit(
				temp_new_dat$B_0_fourth_root
			)[1])
		}
		
		if ( length(na.omit(temp_new_dat$ln_E)) > 0 )
		{
			df_for_BayesTraits$ln_E[counter] <- paste(
				na.omit(temp_new_dat$ln_E), collapse = ','
			)
			df_for_BayesTraits2$ln_E[counter] <- as.character(na.omit(
				temp_new_dat$ln_E
			)[1])
		}
		
		if ( length(na.omit(temp_new_dat$ln_W_op)) > 0 )
		{
			df_for_BayesTraits$ln_W_op[counter] <- paste(
				na.omit(temp_new_dat$ln_W_op), collapse = ','
			)
			df_for_BayesTraits2$ln_W_op[counter] <- as.character(na.omit(
				temp_new_dat$ln_W_op
			)[1])
		}
		
		if ( length(na.omit(temp_new_dat$T_pk_squared)) > 0 )
		{
			df_for_BayesTraits$T_pk_squared[counter] <- paste(
				na.omit(temp_new_dat$T_pk_squared), collapse = ','
			)
			df_for_BayesTraits2$T_pk_squared[counter] <- as.character(na.omit(
				temp_new_dat$T_pk_squared
			)[1])
		}
		
		if ( length(na.omit(temp_new_dat$ln_B_pk)) > 0 )
		{
			df_for_BayesTraits$ln_B_pk[counter] <- paste(
				na.omit(temp_new_dat$ln_B_pk), collapse = ','
			)
			df_for_BayesTraits2$ln_B_pk[counter] <- as.character(na.omit(
				temp_new_dat$ln_B_pk
			)[1])
		}
		
		if ( length(na.omit(temp_new_dat$ln_E_D)) > 0 )
		{
			df_for_BayesTraits$ln_E_D[counter] <- paste(
				na.omit(temp_new_dat$ln_E_D), collapse = ','
			)
			df_for_BayesTraits2$ln_E_D[counter] <- as.character(na.omit(
				temp_new_dat$ln_E_D
			)[1])
		}
		
		counter <- counter + 1
	}
	
	# Write the tree and the other input files required for the analysis.
	tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% unique(df_for_BayesTraits$Species))])
	write.nexus(tree, file = paste(working_dir, 'tree.nex', sep = ''))
	
	###################
	# B_0_fourth_root #
	###################
	
	write.table(
		df_for_BayesTraits2[,c(1,2)], 
		file = paste(
			working_dir, 'dataset_starting_B_0_fourth_root.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	dat_B_0_fourth_root <- df_for_BayesTraits[,c(1,2,3)]
	dat_B_0_fourth_root <- dat_B_0_fourth_root[dat_B_0_fourth_root[,3] != '-',]
	
	write.table(
		dat_B_0_fourth_root, 
		file = paste(
			working_dir, 'dataset_distdata_B_0_fourth_root.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	sink(paste(working_dir, 'command_B_0_fourth_root_1.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_B_0_fourth_root.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile B_0_fourth_root_run1', 
			'Seed 1337', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	sink(paste(working_dir, 'command_B_0_fourth_root_2.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_B_0_fourth_root.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile B_0_fourth_root_run2', 
			'Seed 1604', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	########
	# ln_E #
	########
	
	write.table(
		df_for_BayesTraits2[,c(1,3)], 
		file = paste(
			working_dir, 'dataset_starting_ln_E.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	dat_ln_E <- df_for_BayesTraits[,c(1,2,4)]
	dat_ln_E <- dat_ln_E[dat_ln_E[,3] != '-',]
	
	write.table(
		dat_ln_E, 
		file = paste(
			working_dir, 'dataset_distdata_ln_E.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	sink(paste(working_dir, 'command_ln_E_1.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_ln_E.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile ln_E_run1', 
			'Seed 1337', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	sink(paste(working_dir, 'command_ln_E_2.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_ln_E.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile ln_E_run2', 
			'Seed 1604', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	###########
	# ln_W_op #
	###########
	
	write.table(
		df_for_BayesTraits2[,c(1,4)], 
		file = paste(
			working_dir, 'dataset_starting_ln_W_op.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	dat_ln_W_op <- df_for_BayesTraits[,c(1,2,5)]
	dat_ln_W_op <- dat_ln_W_op[dat_ln_W_op[,3] != '-',]
	
	write.table(
		dat_ln_W_op, 
		file = paste(
			working_dir, 'dataset_distdata_ln_W_op.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	sink(paste(working_dir, 'command_ln_W_op_1.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_ln_W_op.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile ln_W_op_run1', 
			'Seed 1337', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	sink(paste(working_dir, 'command_ln_W_op_2.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_ln_W_op.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile ln_W_op_run2', 
			'Seed 1604', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	################
	# T_pk_squared #
	################
	
	write.table(
		df_for_BayesTraits2[,c(1,5)], 
		file = paste(
			working_dir, 'dataset_starting_T_pk_squared.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	dat_T_pk_squared <- df_for_BayesTraits[,c(1,2,6)]
	dat_T_pk_squared <- dat_T_pk_squared[dat_T_pk_squared[,3] != '-',]
	
	write.table(
		dat_T_pk_squared, 
		file = paste(
			working_dir, 'dataset_distdata_T_pk_squared.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	sink(paste(working_dir, 'command_T_pk_squared_1.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_T_pk_squared.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile T_pk_squared_run1', 
			'Seed 1337', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	sink(paste(working_dir, 'command_T_pk_squared_2.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_T_pk_squared.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile T_pk_squared_run2', 
			'Seed 1604', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	###########
	# ln_B_pk #
	###########
	
	write.table(
		df_for_BayesTraits2[,c(1,6)], 
		file = paste(
			working_dir, 'dataset_starting_ln_B_pk.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	dat_ln_B_pk <- df_for_BayesTraits[,c(1,2,7)]
	dat_ln_B_pk <- dat_ln_B_pk[dat_ln_B_pk[,3] != '-',]
	
	write.table(
		dat_ln_B_pk, 
		file = paste(
			working_dir, 'dataset_distdata_ln_B_pk.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	sink(paste(working_dir, 'command_ln_B_pk_1.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_ln_B_pk.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile ln_B_pk_run1', 
			'Seed 1337', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	sink(paste(working_dir, 'command_ln_B_pk_2.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_ln_B_pk.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile ln_B_pk_run2', 
			'Seed 1604', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	##########
	# ln_E_D #
	##########
	
	write.table(
		df_for_BayesTraits2[,c(1,7)], 
		file = paste(
			working_dir, 'dataset_starting_ln_E_D.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	dat_ln_E_D <- df_for_BayesTraits[,c(1,2,8)]
	dat_ln_E_D <- dat_ln_E_D[dat_ln_E_D[,3] != '-',]
	
	write.table(
		dat_ln_E_D, 
		file = paste(
			working_dir, 'dataset_distdata_ln_E_D.csv', sep = ''
		), 
		sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE
	)
	
	sink(paste(working_dir, 'command_ln_E_D_1.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_ln_E_D.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile ln_E_D_run1', 
			'Seed 1337', 'Run', sep = "\n"
		)
	)
	
	sink()
	
	sink(paste(working_dir, 'command_ln_E_D_2.txt', sep = ''))
	
	cat(
		paste(
			'4', '2', 'Lambda', 'DistData dataset_distdata_ln_E_D.csv', 
			'Burnin 1000000', 'Iterations 10000000', 'LogFile ln_E_D_run2', 
			'Seed 1604', 'Run', sep = "\n"
		)
	)
	
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

# Read the TPC datasets of phytoplankton and prokaryotes.
dat_phytoplankton <- read.csv('../Data/TPC_parameter_estimates_phytoplankton_r_max.csv', stringsAsFactors = FALSE)
dat_prokaryotes <- read.csv('../Data/TPC_parameter_estimates_prokaryotes_r_max.csv', stringsAsFactors = FALSE)

# Put all Cyanobacteria in the phytoplankton dataset.
dat_phytoplankton <- rbind(dat_phytoplankton, dat_prokaryotes[dat_prokaryotes$Species %in% Cyanobacteria,])

# Remove all Cyanobacteria from the prokaryotes dataset.
dat_prokaryotes <- dat_prokaryotes[!(dat_prokaryotes$Species %in% Cyanobacteria),]

#############################
# P H Y T O P L A N K T O N #
#############################

# Prepare the data file.
prepare_dataset(
	dat_phytoplankton, tree,
	working_dir = '../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/'
)

# Run BayesTraits twice for each TPC parameter.
for ( id in 1:2 )
{
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_B_0_fourth_root.csv < ',
			'command_B_0_fourth_root_', id, '.txt', sep = ''
		)
	)
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_ln_E.csv < ',
			'command_ln_E_', id, '.txt', sep = ''
		)
	)
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_ln_W_op.csv < ',
			'command_ln_W_op_', id, '.txt', sep = ''
		)
	)
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_T_pk_squared.csv < ',
			'command_T_pk_squared_', id, '.txt', sep = ''
		)
	)
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_ln_B_pk.csv < ',
			'command_ln_B_pk_', id, '.txt', sep = ''
		)
	)
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_ln_E_D.csv < ',
			'command_ln_E_D_', id, '.txt', sep = ''
		)
	)
}

# Make sure that the 2 chains per TPC parameter have converged.
run1_phyto_B_0_fourth_root <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/B_0_fourth_root_run1.Log.txt',
	skip = 155
)[,2:118]
run2_phyto_B_0_fourth_root <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/B_0_fourth_root_run2.Log.txt',
	skip = 155
)[,2:118]

check_ESS(run1_phyto_B_0_fourth_root, run2_phyto_B_0_fourth_root)
check_convergence(run1_phyto_B_0_fourth_root, run2_phyto_B_0_fourth_root)

run1_phyto_ln_E <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ln_E_run1.Log.txt',
	skip = 155
)[,2:118]
run2_phyto_ln_E <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ln_E_run2.Log.txt',
	skip = 155
)[,2:118]

check_ESS(run1_phyto_ln_E, run2_phyto_ln_E)
check_convergence(run1_phyto_ln_E, run2_phyto_ln_E)

run1_phyto_ln_W_op <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ln_W_op_run1.Log.txt',
	skip = 155
)[,2:118]
run2_phyto_ln_W_op <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ln_W_op_run2.Log.txt',
	skip = 155
)[,2:118]

check_ESS(run1_phyto_ln_W_op, run2_phyto_ln_W_op)
check_convergence(run1_phyto_ln_W_op, run2_phyto_ln_W_op)

run1_phyto_T_pk_squared <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/T_pk_squared_run1.Log.txt',
	skip = 155
)[,2:118]
run2_phyto_T_pk_squared <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/T_pk_squared_run2.Log.txt',
	skip = 155
)[,2:118]

check_ESS(run1_phyto_T_pk_squared, run2_phyto_T_pk_squared)
check_convergence(run1_phyto_T_pk_squared, run2_phyto_T_pk_squared)

run1_phyto_ln_B_pk <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ln_B_pk_run1.Log.txt',
	skip = 155
)[,2:118]
run2_phyto_ln_B_pk <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ln_B_pk_run2.Log.txt',
	skip = 155
)[,2:118]

check_ESS(run1_phyto_ln_B_pk, run2_phyto_ln_B_pk)
check_convergence(run1_phyto_ln_B_pk, run2_phyto_ln_B_pk)

run1_phyto_ln_E_D <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ln_E_D_run1.Log.txt',
	skip = 155
)[,2:118]
run2_phyto_ln_E_D <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/phytoplankton/ln_E_D_run2.Log.txt',
	skip = 155
)[,2:118]

check_ESS(run1_phyto_ln_E_D, run2_phyto_ln_E_D)
check_convergence(run1_phyto_ln_E_D, run2_phyto_ln_E_D)

#########################
# P R O K A R Y O T E S #
#########################

# Prepare the data file.
prepare_dataset(
	dat_prokaryotes, tree,
	working_dir = '../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/'
)

# Run BayesTraits twice for each TPC parameter.
for ( id in 1:2 )
{
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_B_0_fourth_root.csv < ',
			'command_B_0_fourth_root_', id, '.txt', sep = ''
		)
	)
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_ln_E.csv < ',
			'command_ln_E_', id, '.txt', sep = ''
		)
	)
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_ln_W_op.csv < ',
			'command_ln_W_op_', id, '.txt', sep = ''
		)
	)
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_T_pk_squared.csv < ',
			'command_T_pk_squared_', id, '.txt', sep = ''
		)
	)
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_ln_B_pk.csv < ',
			'command_ln_B_pk_', id, '.txt', sep = ''
		)
	)
	system(
		paste(
			'cd ../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ ',
			'&& BayesTraitsV3 tree.nex dataset_starting_ln_E_D.csv < ',
			'command_ln_E_D_', id, '.txt', sep = ''
		)
	)
}

# Make sure that the 2 chains per TPC parameter have converged.
run1_prok_B_0_fourth_root <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/B_0_fourth_root_run1.Log.txt',
	skip = 231
)[,2:194]
run2_prok_B_0_fourth_root <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/B_0_fourth_root_run2.Log.txt',
	skip = 231
)[,2:194]

check_ESS(run1_prok_B_0_fourth_root, run2_prok_B_0_fourth_root)
check_convergence(run1_prok_B_0_fourth_root, run2_prok_B_0_fourth_root)

run1_prok_ln_E <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ln_E_run1.Log.txt',
	skip = 231
)[,2:194]
run2_prok_ln_E <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ln_E_run2.Log.txt',
	skip = 231
)[,2:194]

check_ESS(run1_prok_ln_E, run2_prok_ln_E)
check_convergence(run1_prok_ln_E, run2_prok_ln_E)

run1_prok_ln_W_op <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ln_W_op_run1.Log.txt',
	skip = 231
)[,2:194]
run2_prok_ln_W_op <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ln_W_op_run2.Log.txt',
	skip = 231
)[,2:194]

check_ESS(run1_prok_ln_W_op, run2_prok_ln_W_op)
check_convergence(run1_prok_ln_W_op, run2_prok_ln_W_op)

run1_prok_T_pk_squared <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/T_pk_squared_run1.Log.txt',
	skip = 231
)[,2:194]
run2_prok_T_pk_squared <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/T_pk_squared_run2.Log.txt',
	skip = 231
)[,2:194]

check_ESS(run1_prok_T_pk_squared, run2_prok_T_pk_squared)
check_convergence(run1_prok_T_pk_squared, run2_prok_T_pk_squared)

run1_prok_ln_B_pk <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ln_B_pk_run1.Log.txt',
	skip = 231
)[,2:194]
run2_prok_ln_B_pk <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ln_B_pk_run2.Log.txt',
	skip = 231
)[,2:194]

check_ESS(run1_prok_ln_B_pk, run2_prok_ln_B_pk)
check_convergence(run1_prok_ln_B_pk, run2_prok_ln_B_pk)

run1_prok_ln_E_D <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ln_E_D_run1.Log.txt',
	skip = 231
)[,2:194]
run2_prok_ln_E_D <- read.delim(
	'../Results/phylogenetic_heritabilities_BayesTraits/prokaryotes/ln_E_D_run2.Log.txt',
	skip = 231
)[,2:194]

check_ESS(run1_prok_ln_E_D, run2_prok_ln_E_D)
check_convergence(run1_prok_ln_E_D, run2_prok_ln_E_D)
