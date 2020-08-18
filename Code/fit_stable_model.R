#!/usr/bin/Rscript

# This script is used to estimate the evolutionary rates of thermal 
# sensitivity based on the stable model (Fig 5), and to generate 
# plots of thermal sensitivity evolution through time (Figs 6 and L in 
# the S1 Appendix).

library(ape)
library(coda)
library(geiger)
library(phytools)

#####################
# F U N C T I O N S #
#####################

# This function fits the stable model of trait evolution to ln(E) values. 
# It then generates a plot of the phylogeny, with branches coloured 
# according to the normalised rate of evolution. Finally, it also 
# plots the evolution of trait values from the root until the present 
# time.
fit_stable_model_E <- function(tree, dat, working_dir)
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
	
	# Write the tree and the dataset so that stabletraits can analyse 
	# them.
	write.tree(tree, file = paste(working_dir, "tree.nwk", sep = ''))
	
	write.table(
		dat[,c("Species", "ln_E")], col.names = FALSE, 
		row.names = FALSE, sep = "\t",
		file = paste(working_dir, "dataset.txt", sep = ''), 
		quote = FALSE
	)
	
	# Run stabletraits.
	system(
		paste(
			'stabletraits -r 1337 -o ', working_dir, 'out -d ',
			working_dir, 'dataset.txt -t ', working_dir, 'tree.nwk -c 4 -i 30000000',
			sep = ''
		)
	)
	
	# Make sure that the 4 chains converged and that the Effective Sample 
	# Size for each parameter is high enough.
	chain1 <- read.delim(paste(working_dir, 'out.chain1.log', sep = ''))
	chain2 <- read.delim(paste(working_dir, 'out.chain2.log', sep = ''))
	chain3 <- read.delim(paste(working_dir, 'out.chain3.log', sep = ''))
	chain4 <- read.delim(paste(working_dir, 'out.chain4.log', sep = ''))
	
	chain1 <- chain1[(round(nrow(chain1) * 0.25) + 1):nrow(chain1),2:ncol(chain1)]
	chain2 <- chain2[(round(nrow(chain2) * 0.25) + 1):nrow(chain2),2:ncol(chain2)]
	chain3 <- chain3[(round(nrow(chain3) * 0.25) + 1):nrow(chain3),2:ncol(chain3)]
	chain4 <- chain4[(round(nrow(chain4) * 0.25) + 1):nrow(chain4),2:ncol(chain4)]
	
	if (
		min(effectiveSize(chain1[,1:(ncol(chain1) - 1)])) < 200 ||
		min(effectiveSize(chain2[,1:(ncol(chain2) - 1)])) < 200 ||
		min(effectiveSize(chain3[,1:(ncol(chain3) - 1)])) < 200 ||
		min(effectiveSize(chain4[,1:(ncol(chain4) - 1)])) < 200 ||
		max(gelman.diag(mcmc.list(
					mcmc(chain1[,1:(ncol(chain1) - 1)]), 
					mcmc(chain2[,1:(ncol(chain2) - 1)]), 
					mcmc(chain3[,1:(ncol(chain3) - 1)]), 
					mcmc(chain4[,1:(ncol(chain4) - 1)])
				), multivariate = FALSE
			)$psrf[,1]
		) >= 1.1
	)
	{
		unlink(paste(working_dir, '*', sep = ''))
		stop("Stabletraits hasn't converged! Increase the chain length!")
	}	

	# Process the output of stabletraits.
	system(
		paste(
			'stabletraitssum --path ', working_dir, 'out ',
			'--from 7500000 --to 30000000',
			sep = ''
		)
	)
	
	system(paste('perl reformat_stabletraits_rates.pl ', working_dir, sep = ''))
	
	# Extract the estimated ancestral trait values and the uncertainty 
	# around them.
	all_samples <- rbind(chain1, chain2, chain3, chain4)
	ancestral_states <- exp(colMeans(all_samples[,4:(ncol(all_samples) - 1)]))
	names(ancestral_states) <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)

	HPD_95_high_E <- exp(apply(all_samples[,4:(ncol(all_samples) - 1)], FUN = function(x) HPDinterval(mcmc(x)), MARGIN = 2)[2,])
	names(HPD_95_high_E) <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)
	
	HPD_95_low_E <- exp(apply(all_samples[,4:(ncol(all_samples) - 1)], FUN = function(x) HPDinterval(mcmc(x)), MARGIN = 2)[1,])
	names(HPD_95_low_E) <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)

	sp_estimates <- exp(dat$ln_E)
	names(sp_estimates) <- dat$Species

	# Make a plot of the normalised rate of evolution for each branch 
	# of the phylogeny.
	tree_stabletraits <- read.newick(paste(working_dir, 'tree.nwk', sep = ''))
	dat_rates <- read.delim(paste(working_dir, 'out.brlens_reformatted', sep = ''))
	
	pdf(file = paste(working_dir, 'tree_rates.pdf', sep = ''), width = 2.5, height = 2.5)
	plotBranchbyTrait(
		tree_stabletraits, dat_rates$Rate/max(dat_rates$Rate), 
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

	# Generate a phenogram of trait values.
	pdf(file = paste(working_dir, 'tree.pdf', sep = ''), width = 2.5, height = 2.5)

	par(mar = c(3.25,3.25,0,0), mgp = c(2,0.5,0), cex = 0.6 )

	my_phenogram95(
		rescale(tree, 'depth', 1), HPD_95_high_E, HPD_95_low_E, ancestral_states, x = sp_estimates,
		ylim = c(0,4),
		xlab = "Relative time",
		ylab = bquote(italic(E) ~ "(eV)")
	)
	
	dev.off()
	
	gc()
	
	return()
}

# This function fits the stable model of trait evolution to ln(W_op) values. 
# It then generates a plot of the phylogeny, with branches coloured 
# according to the normalised rate of evolution. Finally, it also 
# plots the evolution of trait values from the root until the present 
# time.
fit_stable_model_W_op <- function(tree, dat, working_dir)
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

	# Write the tree and the dataset so that stabletraits can analyse 
	# them.		
	write.tree(tree, file = paste(working_dir, "tree.nwk", sep = ''))
	
	write.table(
		dat[,c("Species", "ln_W_op")], col.names = FALSE, 
		row.names = FALSE, sep = "\t",
		file = paste(working_dir, "dataset.txt", sep = ''), 
		quote = FALSE
	)

	# Run stabletraits.
	system(
		paste(
			'stabletraits -r 1337 -o ', working_dir, 'out -d ',
			working_dir, 'dataset.txt -t ', working_dir, 'tree.nwk -c 4 -i 30000000',
			sep = ''
		)
	)

	# Make sure that the 4 chains converged and that the Effective Sample 
	# Size for each parameter is high enough.
	chain1 <- read.delim(paste(working_dir, 'out.chain1.log', sep = ''))
	chain2 <- read.delim(paste(working_dir, 'out.chain2.log', sep = ''))
	chain3 <- read.delim(paste(working_dir, 'out.chain3.log', sep = ''))
	chain4 <- read.delim(paste(working_dir, 'out.chain4.log', sep = ''))
	
	chain1 <- chain1[(round(nrow(chain1) * 0.25) + 1):nrow(chain1),2:ncol(chain1)]
	chain2 <- chain2[(round(nrow(chain2) * 0.25) + 1):nrow(chain2),2:ncol(chain2)]
	chain3 <- chain3[(round(nrow(chain3) * 0.25) + 1):nrow(chain3),2:ncol(chain3)]
	chain4 <- chain4[(round(nrow(chain4) * 0.25) + 1):nrow(chain4),2:ncol(chain4)]
	
	if (
		min(effectiveSize(chain1[,1:(ncol(chain1) - 1)])) < 200 ||
		min(effectiveSize(chain2[,1:(ncol(chain2) - 1)])) < 200 ||
		min(effectiveSize(chain3[,1:(ncol(chain3) - 1)])) < 200 ||
		min(effectiveSize(chain4[,1:(ncol(chain4) - 1)])) < 200 ||
		max(gelman.diag(mcmc.list(
					mcmc(chain1[,1:(ncol(chain1) - 1)]), 
					mcmc(chain2[,1:(ncol(chain2) - 1)]), 
					mcmc(chain3[,1:(ncol(chain3) - 1)]), 
					mcmc(chain4[,1:(ncol(chain4) - 1)])
				), multivariate = FALSE
			)$psrf[,1]
		) >= 1.1
	)
	{
		unlink(paste(working_dir, '*', sep = ''))
		stop("Stabletraits hasn't converged! Increase the chain length!")
	}

	# Process the output of stabletraits.
	system(
		paste(
			'stabletraitssum --path ', working_dir, 'out ',
			'--from 7500000 --to 30000000',
			sep = ''
		)
	)
	
	system(paste('perl reformat_stabletraits_rates.pl ', working_dir, sep = ''))

	# Extract the estimated ancestral trait values and the uncertainty 
	# around them.	
	all_samples <- rbind(chain1, chain2, chain3, chain4)
	ancestral_states <- exp(colMeans(all_samples[,4:(ncol(all_samples) - 1)]))
	names(ancestral_states) <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)

	HPD_95_high_W_op <- exp(apply(all_samples[,4:(ncol(all_samples) - 1)], FUN = function(x) HPDinterval(mcmc(x)), MARGIN = 2)[2,])
	names(HPD_95_high_W_op) <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)
	
	HPD_95_low_W_op <- exp(apply(all_samples[,4:(ncol(all_samples) - 1)], FUN = function(x) HPDinterval(mcmc(x)), MARGIN = 2)[1,])
	names(HPD_95_low_W_op) <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)

	sp_estimates <- exp(dat$ln_W_op)
	names(sp_estimates) <- dat$Species

	# Make a plot of the normalised rate of evolution for each branch 
	# of the phylogeny.
	tree_stabletraits <- read.newick(paste(working_dir, 'tree.nwk', sep = ''))
	dat_rates <- read.delim(paste(working_dir, 'out.brlens_reformatted', sep = ''))
	
	pdf(file = paste(working_dir, 'tree_rates.pdf', sep = ''), width = 2.5, height = 2.5)
	plotBranchbyTrait(
		tree_stabletraits, dat_rates$Rate/max(dat_rates$Rate), 
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

	# Generate a phenogram of trait values.
	pdf(file = paste(working_dir, 'tree.pdf', sep = ''), width = 2.5, height = 2.5)
	par(mar = c(3.25,3.25,0,0), mgp = c(2,0.5,0), cex = 0.6 )

	my_phenogram95(
		rescale(tree, 'depth', 1), HPD_95_high_W_op, HPD_95_low_W_op, ancestral_states, x = sp_estimates,
		ylim = c(0,40),
		xlab = "Relative time",
		ylab = bquote(italic(W)[op] ~ "(°C)")
	)
	
	dev.off()
	
	gc()
	
	return()
}

# This function comes from the phytools R package and has been slightly 
# modified for visual purposes. It projects the phylogeny into trait 
# value vs time space.
my_phenogram95 <- function (tree, CI_95_high, CI_95_low, means, ...) 
{
    if (hasArg(x)) 
        x <- list(...)$x
    else stop("no phenotypic data provided")
    if (hasArg(spread.labels)) 
        spread.labels <- list(...)$spread.labels
    else spread.labels <- TRUE
    if (hasArg(link)) 
        link <- list(...)$link
    else link <- 0.05 * max(nodeHeights(tree))
    if (hasArg(offset)) 
        offset <- list(...)$offset
    else offset <- 0
    if (hasArg(hold)) 
        hold <- list(...)$hold
    else hold <- TRUE
    if (hasArg(quiet)) 
        quiet <- list(...)$quiet
    else quiet <- FALSE
    if (hasArg(tlim)) 
        tlim <- list(...)$tlim
    else tlim <- c(0, 25)
    trans <- as.character(floor(seq(tlim[1], tlim[2], length.out = 51)))
    trans[as.numeric(trans) < 10] <- paste("0", trans[as.numeric(trans) < 
        10], sep = "")
    args <- list(...)
    args$tree <- tree
    args$lwd <- 1
    args$link <- 0.05 * max(nodeHeights(tree))
    args$offset <- 0
    if (hold) 
        null <- dev.hold()
    if (!quiet && hold) {
        cat("Computing density traitgram...\n")
        flush.console()
    }
    for (i in 0:50) {
        p <- i/length(trans)
        args$add <- i > 0
        args$ftype <- "off"
        args$spread.labels <- FALSE
        args$link <- if (i == 0) 
            link
        else 0
        args$offset <- if (i == 0) 
            offset
        else offset + link
        args$x <- c(x, (1 - p) * CI_95_high + p * means)
        args$colors <- paste("#dd0000", trans[i + 1], sep = "")
        do.call(phenogram, args)
        args$x <- c(x, (1 - p) * CI_95_low + p * means)
        args$add <- TRUE
        args$link <- 0
        args$offset <- offset + link
        do.call(phenogram, args)
    }
    args$x <- c(x, means)
    args$add <- TRUE
    args$colors <- "#ffd700"
    args$lwd <- 1.2
    args$offset <- offset + link
    args$ftype <- "off"
    args$spread.labels <- FALSE
    do.call(phenogram, args)
    null <- dev.flush()
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
# stable model of trait evolution.
dat_phytoplankton <- rbind(dat_phytoplankton, dat_prokaryotes[dat_prokaryotes$Species %in% Cyanobacteria,])
fit_stable_model_E(
	tree, dat_phytoplankton, '../Results/stable/E/phytoplankton/'
)
fit_stable_model_W_op(
	tree, dat_phytoplankton, '../Results/stable/W_op/phytoplankton/'
)

# Remove all Cyanobacteria from the prokaryotes dataset and fit the 
# stable model of trait evolution.
dat_prokaryotes <- dat_prokaryotes[!(dat_prokaryotes$Species %in% Cyanobacteria),]
fit_stable_model_E(
	tree, dat_prokaryotes, '../Results/stable/E/prokaryotes/'
)
fit_stable_model_W_op(
	tree, dat_prokaryotes, '../Results/stable/W_op/prokaryotes/'
)

# Fit the stable model also to the photosynthesis and respiration datasets.
fit_stable_model_E(
	tree, dat_photosynthesis, '../Results/stable/E/photosynthesis/'
)
fit_stable_model_W_op(
	tree, dat_photosynthesis, '../Results/stable/W_op/photosynthesis/'
)

fit_stable_model_E(
	tree, dat_respiration, '../Results/stable/E/respiration/'
)
fit_stable_model_W_op(
	tree, dat_respiration, '../Results/stable/W_op/respiration/'
)
