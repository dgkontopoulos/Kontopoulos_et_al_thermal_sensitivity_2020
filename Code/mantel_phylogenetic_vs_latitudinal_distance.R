#!/usr/bin/Rscript

# This script performs a Mantel test to estimate the correlation 
# between phylogenetic distance and latitudinal distance for the 
# datasets of phytoplankton and prokaryotes.

library(ade4)
library(ape)

#####################
# F U N C T I O N S #
#####################

# This function runs a Mantel test of phylogenetic vs latitudinal 
# distance. The results are written to a file.
run_mantel_test <- function(dataset, tree, output_file)
{

	# Prepare the dataset.
	dataset$Species <- gsub("\\W", "_", dataset$Species)
	dataset <- dataset[dataset$Species %in% tree$tip.label,]
	dataset <- na.omit(dataset[,c('Species', 'Latitude')])
	dataset <- dataset[!duplicated(dataset),]
	
	# Prepare the latitudinal distance matrix.
	mat_lat <- matrix(nrow = dim(dataset)[1], ncol = dim(dataset)[1])
	diag(mat_lat) <- 0
	
	# Prepare the phylogenetic distance matrix.
	mat_evol <- matrix(nrow = dim(dataset)[1], ncol = dim(dataset)[1])
	diag(mat_evol) <- 0
	
	# Populate the matrices.
	cophenetic_dists <- cophenetic(tree)
	
	for ( i in 1:(nrow(dataset) - 1) )
	{
		for ( j in i:nrow(dataset) )
		{
			mat_evol[i,j] <- cophenetic_dists[
				dataset$Species[i], dataset$Species[j]
			]
			mat_evol[j,i] <- mat_evol[i,j]
	
			mat_lat[i,j] <- sqrt((dataset$Latitude[i] - dataset$Latitude[j])^2)
			mat_lat[j,i] <- mat_lat[i,j]
		}
	}

	# Execute the Mantel test with 9,999 permutations.
	set.seed(1337)
	mantel_res <- mantel.rtest(as.dist(mat_lat), as.dist(mat_evol), nrepet = 9999)
	
	# Write the results of the test to the output file.
	sink(file = output_file)
	print(mantel_res)
	sink(NULL)
	
	return(0)
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

# Put all Cyanobacteria in the phytoplankton dataset and run the Mantel test.
dat_phytoplankton <- rbind(dat_phytoplankton, dat_prokaryotes[dat_prokaryotes$Species %in% Cyanobacteria,])
run_mantel_test(dat_phytoplankton, tree, '../Results/Mantel_phylogenetic_vs_latitude_distance_phytoplankton.txt')

# Remove all Cyanobacteria from the prokaryotes dataset and run the Mantel test.
dat_prokaryotes <- dat_prokaryotes[!(dat_prokaryotes$Species %in% Cyanobacteria),]
run_mantel_test(dat_prokaryotes, tree, '../Results/Mantel_phylogenetic_vs_latitude_distance_prokaryotes.txt')
