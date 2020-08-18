#!/usr/bin/Rscript

# This script fits the free model of trait evolution to ln(E) and 
# ln(W_op) estimates of phytoplankton and prokaryotes.
#
# The output is a PDF plot of the tree, with branches coloured 
# according to the normalised evolutionary rate (Fig G in the S1 
# Appendix).

library(ape)
library(motmot.2.0)
library(phytools)

#####################
# F U N C T I O N S #
#####################

# This function fits the free model to the ln(E) values in a dataset.
run_free_rates_E <- function(tree, dat, output_file)
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
	
	# Fit the free model.
	dat_for_free <- dat$ln_E
	names(dat_for_free) <- dat$Species
	
	fit_all_rates <- transformPhylo.ML(
		as.matrix(dat_for_free), tree, model = 'free'
	)
	
	# Plot the phylogeny with branches coloured according to the 
	# normalised rate of evolution (normalised sigma-squared).
	pdf(file = output_file, width = 2.5, height = 2.5)
	plotBranchbyTrait(
		tree, fit_all_rates$Rates/max(fit_all_rates$Rates), 
		type = 'fan', cex = 0.001, edge.width = 2, legend = FALSE,
		palette = colorRampPalette(
			c(
				colorRampPalette(c('#ffeda0', '#fed976'))(43), 
				colorRampPalette(c('#fed976', '#feb24c'))(67),
				colorRampPalette(c('#feb24c', '#fd8d3c'))(100),
				colorRampPalette(c('#fd8d3c', '#fc4e2a'))(117),
				colorRampPalette(c('#fc4e2a', '#e31a1c'))(143),
				colorRampPalette(c('#e31a1c', '#bd0026'))(271),
				colorRampPalette(c('#bd0026', '#800026'))(259)
			)
		)
	)
	dev.off()
	
	return()
}

# This function fits the free model to the ln(W_op) values in a dataset.
run_free_rates_W_op <- function(tree, dat, output_file)
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

	# Fit the free model.	
	dat_for_free <- dat$ln_W_op
	names(dat_for_free) <- dat$Species
	
	fit_all_rates <- transformPhylo.ML(
		as.matrix(dat_for_free), tree, model = 'free'
	)

	# Plot the phylogeny with branches coloured according to the 
	# normalised rate of evolution (normalised sigma-squared).	
	pdf(file = output_file, width = 2.5, height = 2.5)
	plotBranchbyTrait(
		tree, fit_all_rates$Rates/max(fit_all_rates$Rates), 
		type = 'fan', cex = 0.001, edge.width = 2, legend = FALSE,
		palette = colorRampPalette(
			c(
				colorRampPalette(c('#ffeda0', '#fed976'))(43), 
				colorRampPalette(c('#fed976', '#feb24c'))(67),
				colorRampPalette(c('#feb24c', '#fd8d3c'))(100),
				colorRampPalette(c('#fd8d3c', '#fc4e2a'))(117),
				colorRampPalette(c('#fc4e2a', '#e31a1c'))(143),
				colorRampPalette(c('#e31a1c', '#bd0026'))(271),
				colorRampPalette(c('#bd0026', '#800026'))(259)
			)
		)
	)
	dev.off()
	
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

# Put all Cyanobacteria in the phytoplankton dataset and fit the free model for E and W_op.
dat_phytoplankton <- rbind(dat_phytoplankton, dat_prokaryotes[dat_prokaryotes$Species %in% Cyanobacteria,])
run_free_rates_E(
	tree, dat_phytoplankton, 
	'../Results/free_rates_E_phytoplankton.pdf'
)
run_free_rates_W_op(
	tree, dat_phytoplankton, 
	'../Results/free_rates_W_op_phytoplankton.pdf'
)

# Remove all Cyanobacteria from the prokaryotes dataset and fit the free model for E and W_op.
dat_prokaryotes <- dat_prokaryotes[!(dat_prokaryotes$Species %in% Cyanobacteria),]
run_free_rates_E(
	tree, dat_prokaryotes, 
	'../Results/free_rates_E_prokaryotes.pdf'
)
run_free_rates_W_op(
	tree, dat_prokaryotes, 
	'../Results/free_rates_W_op_prokaryotes.pdf'
)
