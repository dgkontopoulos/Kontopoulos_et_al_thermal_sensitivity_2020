#!/usr/bin/Rscript

# This script performs disparity-through-time analyses for ln(E) and 
# ln(W_op) estimates of phytoplankton (with or without Cyanobacteria) 
# and prokaryotes.

library(ape)
library(geiger)
library(spptest)

#####################
# F U N C T I O N S #
#####################

# This function performs a disparity-through-time analysis for ln(E) 
# using the rank envelope test approach of David J. Murrell.
run_dtt_envelope_E <- function(tree, dat, output_file, max_disp = 2.6)
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
	dat_for_dtt <- dat$ln_E
	names(dat_for_dtt) <- dat$Species
	
	# The rest of this function is mostly based on D. J. Murrell's code 
	# for performing a disparity-through-time analysis using a rank 
	# envelope test (see https://github.com/djmurrell/DTT-Envelope-code/blob/master/empirical.tests.R).
	
	source("https://raw.githubusercontent.com/mwpennell/geiger-v2/master/R/disparity.R")
	source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/dtt1.R")
	source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/getMDI1.R")
	source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/rank_dtt.R")

	# Run the analysis using 10,000 simulations of Brownian motion.
	set.seed(1337)
	d1 <- dtt1(rescale(tree, 'depth', 1), dat_for_dtt, nsim = 10000, plot = T, Ylim = c(0,max_disp), calculateMDIp = TRUE)
	
	ylim<-par("yaxp")
	r1<-rank_env_dtt(d1, Plot=F)

	# Print the P-value range.
	print(r1$p_interval)
	
	# Plot the result.
	pdf(output_file, width = 2.5, height = 2.5)

	par(
		mai=c(0.5,0.45,0.01,0.01), cex.main = 0.75, cex.axis = 0.6,
        cex.lab = 0.75, cex.sub = 0.1, mgp = c(1.5, 0.5, 0)
	)
	
	plot(
		c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), 
		type="n", xlab="Relative time", frame.plot=F, 
		ylab=bquote("Mean subclade disparity in ln(" * italic(E) * ")"), main= ''
	)
	x<-r1$r
	y1<-r1$upper
	y2<-r1$lower
	polygon(c(x, rev(x)), c(y1, rev(y2)), col="grey60", border=NA)
	lines(x[1:(length(x) - 1)], r1$data_curve[1:(length(x) - 1)], lwd=2)
	lines(x, r1$central_curve, lty=2)
	
	dev.off()
	
	return()
}

# This function performs a disparity-through-time analysis for ln(W_op) 
# using the rank envelope test approach of David J. Murrell.
run_dtt_envelope_W_op <- function(tree, dat, output_file, max_disp = 3)
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
	dat_for_dtt <- dat$ln_W_op
	names(dat_for_dtt) <- dat$Species
	
	# The rest of this function is mostly based on D. J. Murrell's code 
	# for performing a disparity-through-time analysis using a rank 
	# envelope test (see https://github.com/djmurrell/DTT-Envelope-code/blob/master/empirical.tests.R).
	
	source("https://raw.githubusercontent.com/mwpennell/geiger-v2/master/R/disparity.R")
	source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/dtt1.R")
	source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/getMDI1.R")
	source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/rank_dtt.R")

	# Run the analysis using 10,000 simulations of Brownian motion.
	set.seed(1337)
	d1 <- dtt1(rescale(tree, 'depth', 1), dat_for_dtt, nsim = 10000, plot = T, Ylim = c(0,max_disp), calculateMDIp = TRUE)
	
	ylim<-par("yaxp")
	r1<-rank_env_dtt(d1, Plot=F)

	# Print the P-value range.
	print(r1$p_interval)
	
	# Plot the result.
	pdf(output_file, width = 2.5, height = 2.5)

	par(
		mai=c(0.5,0.45,0.01,0.01), cex.main = 0.75, cex.axis = 0.6,
        cex.lab = 0.75, cex.sub = 0.1, mgp = c(1.5, 0.5, 0)
	)
	
	plot(
		c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), 
		type="n", xlab="Relative time", frame.plot=F, 
		ylab=bquote("Mean subclade disparity in ln(" * italic(E) * ")"), main= ''
	)
	x<-r1$r
	y1<-r1$upper
	y2<-r1$lower
	polygon(c(x, rev(x)), c(y1, rev(y2)), col="grey60", border=NA)
	lines(x[1:(length(x) - 1)], r1$data_curve[1:(length(x) - 1)], lwd=2)
	lines(x, r1$central_curve, lty=2)
	
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

# Remove Cyanobacteria from the phytoplankton dataset and perform a 
# disparity-through-time analysis for ln(E) and ln(W_op).
dat_eukaryotic_phytoplankton <- dat_phytoplankton[!(dat_phytoplankton$Species %in% Cyanobacteria),]
run_dtt_envelope_E(
	tree, dat_eukaryotic_phytoplankton,
	'../Results/dtt_E_eukaryotic_phytoplankton.pdf',
	max_disp = 2.5
)
run_dtt_envelope_W_op(
	tree, dat_eukaryotic_phytoplankton,
	'../Results/dtt_W_op_eukaryotic_phytoplankton.pdf',
	max_disp = 2
)

# Put all Cyanobacteria in the phytoplankton dataset and perform a 
# disparity-through-time analysis for ln(E) and ln(W_op).
dat_phytoplankton <- rbind(dat_phytoplankton, dat_prokaryotes[dat_prokaryotes$Species %in% Cyanobacteria,])
run_dtt_envelope_E(
	tree, dat_phytoplankton,
	'../Results/dtt_E_phytoplankton.pdf',
	max_disp = 2.5
)
run_dtt_envelope_W_op(
	tree, dat_phytoplankton,
	'../Results/dtt_W_op_phytoplankton.pdf',
	max_disp = 3.4
)

# Remove all Cyanobacteria from the prokaryotes dataset and perform a 
# disparity-through-time analysis for ln(E) and ln(W_op).
dat_prokaryotes <- dat_prokaryotes[!(dat_prokaryotes$Species %in% Cyanobacteria),]
run_dtt_envelope_E(
	tree, dat_prokaryotes,
	'../Results/dtt_E_prokaryotes.pdf',
	max_disp = 2
)
run_dtt_envelope_W_op(
	tree, dat_prokaryotes,
	'../Results/dtt_W_op_prokaryotes.pdf',
	max_disp = 3
)
