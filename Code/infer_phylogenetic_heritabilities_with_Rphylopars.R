#!/usr/bin/Rscript

# This script estimates the phylogenetic heritability (Pagel's lambda) 
# separately for each TPC parameter of phytoplankton and prokaryotes, 
# using the Rphylopars R package.
#
# The results are written to an output file.

library(Rphylopars)

#####################
# F U N C T I O N S #
#####################

# This function prepares TPC datasets (e.g., removes spaces in species 
# names, removes estimates from bad fits of the Sharpe-Schoolfield 
# model etc.) for the estimation of phylogenetic heritabilities 
# (Pagel's lambda). 
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

#############################
# P H Y T O P L A N K T O N #
#############################

# Put all Cyanobacteria in the phytoplankton dataset and estimate 
# phylogenetic heritabilities.
dat_phytoplankton <- rbind(dat_phytoplankton, dat_prokaryotes[dat_prokaryotes$Species %in% Cyanobacteria,])
dat_phytoplankton_cleaned <- prepare_dataset(dat_phytoplankton, tree)

dat_for_phylopars_phytoplankton <- dat_phytoplankton_cleaned$dat[,c(
		'Species',	'B_0_fourth_root',
		'ln_E',				'ln_W_op',
		'T_pk_squared',		'ln_B_pk',
		'ln_E_D'
	)
]

names(dat_for_phylopars_phytoplankton)[1] <- 'species'

phytoplankton_lambda_B_0_fourth_root <- phylopars(
	dat_for_phylopars_phytoplankton[,c(1,2)], dat_phytoplankton_cleaned$tree,
	model = 'lambda'
)$model$lambda

phytoplankton_lambda_ln_E <- phylopars(
	dat_for_phylopars_phytoplankton[,c(1,3)], dat_phytoplankton_cleaned$tree,
	model = 'lambda'
)$model$lambda

phytoplankton_lambda_ln_W_op <- phylopars(
	dat_for_phylopars_phytoplankton[,c(1,4)], dat_phytoplankton_cleaned$tree,
	model = 'lambda'
)$model$lambda

phytoplankton_lambda_T_pk_squared <- phylopars(
	dat_for_phylopars_phytoplankton[,c(1,5)], dat_phytoplankton_cleaned$tree,
	model = 'lambda'
)$model$lambda

phytoplankton_lambda_ln_B_pk <- phylopars(
	dat_for_phylopars_phytoplankton[,c(1,6)], dat_phytoplankton_cleaned$tree,
	model = 'lambda'
)$model$lambda

phytoplankton_lambda_ln_E_D <- phylopars(
	dat_for_phylopars_phytoplankton[,c(1,7)], dat_phytoplankton_cleaned$tree,
	model = 'lambda'
)$model$lambda

#########################
# P R O K A R Y O T E S #
#########################

# Remove all Cyanobacteria from the prokaryotes dataset and estimate 
# phylogenetic heritabilities.
dat_prokaryotes <- dat_prokaryotes[!(dat_prokaryotes$Species %in% Cyanobacteria),]
dat_prokaryotes_cleaned <- prepare_dataset(dat_prokaryotes, tree)

dat_for_phylopars_prokaryotes <- dat_prokaryotes_cleaned$dat[,c(
		'Species',	'B_0_fourth_root',
		'ln_E',				'ln_W_op',
		'T_pk_squared',		'ln_B_pk',
		'ln_E_D'
	)
]

names(dat_for_phylopars_prokaryotes)[1] <- 'species'

prokaryotes_lambda_B_0_fourth_root <- phylopars(
	dat_for_phylopars_prokaryotes[,c(1,2)], dat_prokaryotes_cleaned$tree,
	model = 'lambda'
)$model$lambda

prokaryotes_lambda_ln_E <- phylopars(
	dat_for_phylopars_prokaryotes[,c(1,3)], dat_prokaryotes_cleaned$tree,
	model = 'lambda'
)$model$lambda

prokaryotes_lambda_ln_W_op <- phylopars(
	dat_for_phylopars_prokaryotes[,c(1,4)], dat_prokaryotes_cleaned$tree,
	model = 'lambda'
)$model$lambda

prokaryotes_lambda_T_pk_squared <- phylopars(
	dat_for_phylopars_prokaryotes[,c(1,5)], dat_prokaryotes_cleaned$tree,
	model = 'lambda'
)$model$lambda

prokaryotes_lambda_ln_B_pk <- phylopars(
	dat_for_phylopars_prokaryotes[,c(1,6)], dat_prokaryotes_cleaned$tree,
	model = 'lambda'
)$model$lambda

prokaryotes_lambda_ln_E_D <- phylopars(
	dat_for_phylopars_prokaryotes[,c(1,7)], dat_prokaryotes_cleaned$tree,
	model = 'lambda'
)$model$lambda

# Write the resulting lambda estimates to a file.
lambda_df <- data.frame(
	Parameter = rep(
		c(
			'B_0_fourth_root',	'ln_E', 
			'ln_W_op',		'T_pk_squared',
			'ln_B_pk',		'ln_E_D'
		), 2
	),
	Group = c(
		rep('Phytoplankton', 6), rep('Prokaryotes', 6)
	),
	Value = c(
		phytoplankton_lambda_B_0_fourth_root,	phytoplankton_lambda_ln_E,
		phytoplankton_lambda_ln_W_op,			phytoplankton_lambda_T_pk_squared,
		phytoplankton_lambda_ln_B_pk,			phytoplankton_lambda_ln_E_D,	
		prokaryotes_lambda_B_0_fourth_root,		prokaryotes_lambda_ln_E,
		prokaryotes_lambda_ln_W_op,				prokaryotes_lambda_T_pk_squared,
		prokaryotes_lambda_ln_B_pk,				prokaryotes_lambda_ln_E_D
	)
)

write.csv(lambda_df, file = '../Results/phylogenetic_heritabilities_Rphylopars.csv', row.names = FALSE)
