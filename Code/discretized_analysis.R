#!/usr/bin/Rscript

# This script performs the analysis of transitions in the discretized 
# parameter space of E, W_op, and T_pk (Fig I in the S1 Appendix).

library(ape)
library(BAMMtools)
library(phytools)

#####################
# F U N C T I O N S #
#####################

# This function splits the E estimates into 4 bins and fits the Mk model 
# to them. The range for each bin is determined using the Jenks natural 
# breaks clustering algorithm.
fit_Mk_model_E <- function(tree, dat)
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
	
	# Run the Jenks natural breaks algorithm, bin the dataset, and 
	# fit the Mk model.
	dat <- dat[as.character(ids_to_keep),]
	dat_E <- exp(dat$ln_E)
	names(dat_E) <- dat$Species

	my_breaks <- getJenksBreaks(dat_E, 5)
	print(my_breaks)
	
	dat_E_discrete <- dat_E
	dat_E_discrete[dat_E < my_breaks[2]] <- 'A' 
	dat_E_discrete[dat_E >= my_breaks[2] & dat_E < my_breaks[3]] <- 'B'
	dat_E_discrete[dat_E >= my_breaks[3] & dat_E < my_breaks[4]] <- 'C'
	dat_E_discrete[dat_E >= my_breaks[4]] <- 'D'
	
	names(dat_E_discrete) <- dat$Species
	fitted_model <- fitMk(tree, dat_E_discrete, model = 'ARD')
	print(fitted_model)
		
	return(fitted_model)
}

# This function splits the T_pk estimates into 4 bins and fits the Mk model 
# to them. The range for each bin is determined using the Jenks natural 
# breaks clustering algorithm.
fit_Mk_model_T_pk <- function(tree, dat)
{

	# Prepare the dataset (e.g., remove spaces in species names, 
	# remove estimates from bad fits of the Sharpe-Schoolfield model etc.).
	dat$Species <- gsub("\\W", "_", dat$Species)
	dat$Species <- gsub("_+$", "", dat$Species)
	
	dat <- dat[dat$R_squared >= 0.5 & dat$Data_points_rise >= 2 & dat$Data_points_fall >= 2 & dat$T_pk_squared_error_variance != "NaN",]
	
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
	
	# Run the Jenks natural breaks algorithm, bin the dataset, and 
	# fit the Mk model.
	dat <- dat[as.character(ids_to_keep),]
	dat_T_pk <- sqrt(dat$T_pk_squared)
	names(dat_T_pk) <- dat$Species
		
	my_breaks <- getJenksBreaks(dat_T_pk, 5)
	print(my_breaks)
	
	dat_T_pk_discrete <- dat_T_pk
	dat_T_pk_discrete[dat_T_pk < my_breaks[2]] <- 'A' 
	dat_T_pk_discrete[dat_T_pk >= my_breaks[2] & dat_T_pk < my_breaks[3]] <- 'B'
	dat_T_pk_discrete[dat_T_pk >= my_breaks[3] & dat_T_pk < my_breaks[4]] <- 'C'
	dat_T_pk_discrete[dat_T_pk >= my_breaks[4]] <- 'D'
	
	names(dat_T_pk_discrete) <- dat$Species
	fitted_model <- fitMk(tree, dat_T_pk_discrete, model = 'ARD')
	print(fitted_model)
		
	return(fitted_model)
}

# This function splits the W_op estimates into 4 bins and fits the Mk model 
# to them. The range for each bin is determined using the Jenks natural 
# breaks clustering algorithm.
fit_Mk_model_W_op <- function(tree, dat)
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
	
	# Run the Jenks natural breaks algorithm, bin the dataset, and 
	# fit the Mk model.
	dat <- dat[as.character(ids_to_keep),]
	dat_W_op <- exp(dat$ln_W_op)
	names(dat_W_op) <- dat$Species
	
	my_breaks <- getJenksBreaks(dat_W_op, 5)
	print(my_breaks)
	
	dat_W_op_discrete <- dat_W_op
	dat_W_op_discrete[dat_W_op < my_breaks[2]] <- 'A' 
	dat_W_op_discrete[dat_W_op >= my_breaks[2] & dat_W_op < my_breaks[3]] <- 'B'
	dat_W_op_discrete[dat_W_op >= my_breaks[3] & dat_W_op < my_breaks[4]] <- 'C'
	dat_W_op_discrete[dat_W_op >= my_breaks[4]] <- 'D'
	
	names(dat_W_op_discrete) <- dat$Species
	fitted_model <- fitMk(tree, dat_W_op_discrete, model = 'ARD')
	print(fitted_model)
			
	return(fitted_model)
}


####################
# M A I N  C O D E #
####################

# Read the time-calibrated phylogeny.
tree <- read.nexus('../Data/final_calibrated_phylogeny.nex')
tree$node.label <- NULL

# Read the TPC datasets of phytoplankton and prokaryotes.
dat_phytoplankton <- read.csv('../Data/TPC_parameter_estimates_phytoplankton_r_max.csv', stringsAsFactors = FALSE)
dat_prokaryotes <- read.csv('../Data/TPC_parameter_estimates_prokaryotes_r_max.csv', stringsAsFactors = FALSE)

# Merge the two datasets and fit the Mk model to E, W_op, and T_pk.
dat_all <- rbind(dat_phytoplankton, dat_prokaryotes)

E_model <- fit_Mk_model_E(tree, dat_all)
W_op_model <- fit_Mk_model_W_op(tree, dat_all)
T_pk_model <- fit_Mk_model_T_pk(tree, dat_all)
