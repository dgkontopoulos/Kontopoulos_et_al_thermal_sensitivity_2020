#!/usr/bin/Rscript

# This script fits the Lévy model of trait evolution to ln(E) and 
# ln(W_op) estimates of phytoplankton and prokaryotes.
#
# The output is a PDF plot of the tree, with branches coloured 
# according to the posterior probability for the occurrence of a jump 
# (Fig H in the S1 Appendix).

library(phytools)

#####################
# F U N C T I O N S #
#####################

# This function fits the Lévy model of trait evolution to ln(E) values. 
# It then generates a plot of the phylogeny, with branches coloured 
# according to the posterior probability for a jump in trait values.
fit_levy_model_E <- function(tree, dat, working_dir)
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
	
	# Create a directory (if not already available) to store the 
	# output files.
	dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
	
	# Write the tree and the dataset so that levolution can analyse 
	# them.
	write.tree(tree, file = paste(working_dir, "tree.nwk", sep = ''))
	
	write.table(
		dat[,c("Species", "ln_E")], col.names = FALSE, 
		row.names = FALSE, sep = "\t",
		file = paste(working_dir, "dataset.txt", sep = ''), 
		quote = FALSE
	)
	
	# Run levolution.
	system(
		paste(
			'cd ', working_dir, ' && levolution task=infer ', 
			'tree=tree.nwk traits=dataset.txt logAlphaStart=0.5 ', 
			'logAlphaStep=0.5 numOptimizations=5 verbose inferJumps ',
			'maxIterations=2000', sep = ''
		)
	)
	
	# Read the results and generate a plot.
	tree_levy <- read.newick(paste(working_dir, 'tree.nwk', sep = ''))
	dat_levy <- read.tree(paste(working_dir, 'levolution.post', sep = ''))$edge.length
	
	pdf(
		file = paste(working_dir, 'tree_levy.pdf', sep = ''), 
		width = 2.5, height = 2.5
	)
	plotBranchbyTrait(
		tree_levy, dat_levy, 
		type = 'fan', cex = 0.001, edge.width = 2, legend = FALSE,
		palette = colorRampPalette(
			c(
				colorRampPalette(
					c(
						'#000000', '#252525', '#525252', '#737373', 
						'#969696', '#bdbdbd', '#d9d9d9', '#f0f0f0'
					)
				)(750), 
				colorRampPalette(
					c(
						'#ffeda0', '#fed976', '#feb24c', '#fd8d3c', 
						'#fc4e2a', '#e31a1c'
					)
				)(250)
			)
		)
	)
	dev.off()
		
	return()
}

# This function fits the Lévy model of trait evolution to ln(W_op) values. 
# It then generates a plot of the phylogeny, with branches coloured 
# according to the posterior probability for a jump in trait values.
fit_levy_model_W_op <- function(tree, dat, working_dir)
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

	# Create a directory (if not already available) to store the 
	# output files.
	dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)

	# Write the tree and the dataset so that levolution can analyse 
	# them.		
	write.tree(tree, file = paste(working_dir, "tree.nwk", sep = ''))
	
	write.table(
		dat[,c("Species", "ln_W_op")], col.names = FALSE, 
		row.names = FALSE, sep = "\t",
		file = paste(working_dir, "dataset.txt", sep = ''), 
		quote = FALSE
	)

	# Run levolution.
	system(
		paste(
			'cd ', working_dir, ' && levolution task=infer ', 
			'tree=tree.nwk traits=dataset.txt logAlphaStart=0.5 ', 
			'logAlphaStep=0.5 numOptimizations=5 verbose inferJumps ',
			'maxIterations=2000', sep = ''
		)
	)
	
	# Read the results and generate a plot.
	tree_levy <- read.newick(paste(working_dir, 'tree.nwk', sep = ''))
	dat_levy <- read.tree(paste(working_dir, 'levolution.post', sep = ''))$edge.length
	
	pdf(
		file = paste(working_dir, 'tree_levy.pdf', sep = ''), 
		width = 2.5, height = 2.5
	)
	plotBranchbyTrait(
		tree_levy, dat_levy, 
		type = 'fan', cex = 0.001, edge.width = 2, legend = FALSE,
		palette = colorRampPalette(
			c(
				colorRampPalette(
					c(
						'#000000', '#252525', '#525252', '#737373', 
						'#969696', '#bdbdbd', '#d9d9d9', '#f0f0f0'
					)
				)(750), 
				colorRampPalette(
					c(
						'#ffeda0', '#fed976', '#feb24c', '#fd8d3c', 
						'#fc4e2a', '#e31a1c'
					)
				)(250)
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

# Put all Cyanobacteria in the phytoplankton dataset and fit the Lévy model for E and W_op.
dat_phytoplankton <- rbind(dat_phytoplankton, dat_prokaryotes[dat_prokaryotes$Species %in% Cyanobacteria,])
fit_levy_model_E(
	tree, dat_phytoplankton, '../Results/levy/E/phytoplankton/'
)
fit_levy_model_W_op(
	tree, dat_phytoplankton, '../Results/levy/W_op/phytoplankton/'
)

# Remove all Cyanobacteria from the prokaryotes dataset and fit the Lévy model for E and W_op.
dat_prokaryotes <- dat_prokaryotes[!(dat_prokaryotes$Species %in% Cyanobacteria),]
fit_levy_model_E(
	tree, dat_prokaryotes, '../Results/levy/E/prokaryotes/'
)
fit_levy_model_W_op(
	tree, dat_prokaryotes, '../Results/levy/W_op/prokaryotes/'
)
