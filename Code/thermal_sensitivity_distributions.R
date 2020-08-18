#!/usr/bin/Rscript

# This script generates plots of the distributions of E and W_op for 
# phytoplankton and prokaryotes (Figs 4 and C in the S1 Appendix).

library(ape)
library(cowplot)

#####################
# F U N C T I O N S #
#####################

# This function prepares the E datasets (e.g., removes estimates from 
# bad fits of the Sharpe-Schoolfield model) for plotting.
prepare_E_vals <- function(dat, tree)
{
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
	
	return(dat[,c("ln_E", "ln_E_error_variance", "Species")])
}

# This function prepares the W_op datasets (e.g., removes estimates from 
# bad fits of the Sharpe-Schoolfield model) for plotting.
prepare_W_op_vals <- function(dat, tree)
{
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
	
	return(dat[,c("ln_W_op", "ln_W_op_error_variance", "Species")])
}

####################
# M A I N  C O D E #
####################

# Read the time-calibrated phylogeny.
tree <- read.nexus('../Data/final_calibrated_phylogeny.nex')
tree$node.label <- NULL

# Read the TPC datasets of phytoplankton and prokaryotes and clean them.
dat_phytoplankton <- read.csv('../Data/TPC_parameter_estimates_phytoplankton_r_max.csv', stringsAsFactors = FALSE)
dat_prokaryotes <- read.csv('../Data/TPC_parameter_estimates_prokaryotes_r_max.csv', stringsAsFactors = FALSE)

dat_E_phytoplankton_cleaned <- prepare_E_vals(dat_phytoplankton, tree)
dat_E_prokaryotes_cleaned <- prepare_E_vals(dat_prokaryotes, tree)

dat_W_op_phytoplankton_cleaned <- prepare_W_op_vals(dat_phytoplankton, tree)
dat_W_op_prokaryotes_cleaned <- prepare_W_op_vals(dat_prokaryotes, tree)

# Merge the datasets from the two groups.
all_dat_E <- rbind(dat_E_phytoplankton_cleaned, dat_E_prokaryotes_cleaned)
all_dat_W_op <- rbind(dat_W_op_phytoplankton_cleaned, dat_W_op_prokaryotes_cleaned)

################################################
# Species in the largest phyla of our datasets #
################################################

Bacillariophyta <- c(
	'Amphiprora_paludosa', 'Asterionella_formosa',
	'Asterionellopsis_glacialis', 'Coscinodiscus_concinnus',
	'Cyclotella_cryptica', 'Cyclotella_meneghiniana',
	'Ditylum_brightwellii', 'Eucampia_zodiacus',
	'Fragilaria_crotonensis', 'Odontella_aurita',
	'Odontella_sinensis', 'Phaeodactylum_tricornutum',
	'Pseudo_nitzschia_multiseries', 'Pseudo_nitzschia_seriata',
	'Rhizosolenia_setigera', 'Skeletonema_ardens', 
	'Skeletonema_costatum', 'Skeletonema_marinoi',
	'Skeletonema_pseudocostatum', 'Skeletonema_tropicum',
	'Stephanodiscus_hantzschii', 'Thalassionema_nitzschioides',
	'Thalassiosira_eccentrica', 'Thalassiosira_nordenskioeldii',
	'Thalassiosira_pseudonana', 'Thalassiosira_rotula',
	'Thalassiosira_weissflogii' 
)

Cyanobacteria <- c(
	'Mastigocladus_laminosus', 'Anabaena_ucrainica',
	'Aphanizomenon_flosaquae', 'Aphanizomenon_gracile',
	'Aphanizomenon_ovalisporum', 'Cylindrospermopsis_raciborskii',
	'Limnothrix_redekei', 'Microcystis_aeruginosa',
	'Planktothrix_agardhii', 'Prochlorococcus_marinus',
	'Sphaerospermopsis_aphanizomenoides', 'Spirulina_platensis',
	'Synechococcus_elongatus', 'Synechococcus_lividus',
	'Trichodesmium_erythraeum', 'Tychonema_bourrellyi'
)


Dinophyta <- c(
	'Alexandrium_ostenfeldii', 'Alexandrium_tamarense', 'Amphidinium_klebsii',
	'Ceratium_furca', 'Ceratium_fusus', 'Cochlodinium_polykrikoides', 'Coolia_monotis',
	'Gambierdiscus_toxicus', 'Gymnodinium_breve', 'Gymnodinium_catenatum',
	'Gymnodinium_mikimotoi', 'Gymnodinium_veneficum', 'Gyrodinium_instriatum',
	'Heterocapsa_triquetra', 'Prorocentrum_concavum', 'Prorocentrum_dentatum',
	'Prorocentrum_lima', 'Prorocentrum_micans', 'Prorocentrum_minimum',
	'Pyrodinium_bahamense', 'Scrippsiella_trochoidea'
)

Euryarchaeota <- c(
	'Archaeoglobus_veneficus', 'Geoglobus_ahangari', 'Haloarcula_vallismortis',
	'Halobacterium_salinarum', 'Halobaculum_gomorrense', 'Halococcus_morrhuae',
	'Haloferax_volcanii', 'Halogeometricum_borinquense', 'Halorubrum_saccharovorum',
	'Haloterrigena_turkmenica', 'Methanobacterium_subterraneum',
	'Methanococcus_jannaschii', 'Methanococcus_thermolithotrophicus',
	'Methanococcus_voltae', 'Methanoculleus_submarinus', 
	'Methanogenium_frigidum', 'Methanohalophilus_portucalensis',
	'Methanopyrus_kandleri', 'Methanothermococcus_okinawensis',
	'Natrialba_asiatica', 'Natrinema_pellirubrum', 'Natronobacterium_gregoryi',
	'Natronococcus_occultus', 'Natronomonas_pharaonis', 'Natronorubrum_bangense',
	'Palaeococcus_helgesonii', 'Pyrococcus_furiosus', 'Pyrococcus_glycovorans',
	'Thermococcus_alcaliphilus', 'Thermococcus_chitonophagus',
	'Thermococcus_hydrothermalis', 'Thermococcus_siculi', 
	'Thermoplasma_acidophila'
)

Firmicutes <- c(
	'Acetobacterium_carbinolicum', 'Acetobacterium_paludosum',
	'Acetogenium_kivui', 'Aeribacillus_pallidus',
	'Alkaliphilus_transvaalensis', 'Bacillus_acidocaldarius',
	'Bacillus_caldotenax', 'Bacillus_cereus', 'Bacillus_infernus',
	'Bacillus_megaterium', 'Bacillus_subtilis', 'Brochothrix_thermosphacta',
	'Caloramator_indicus', 'Caloranaerobacter_azorensis',
	'Clostridium_autoethanogenum', 'Clostridium_fervidus', 
	'Clostridium_paradoxum', 'Clostridium_perfringens', 'Clostridium_thermoalcaliphilum',
	'Clostridium_thermohydrosulfuricum', 'Clostridium_thermosuccinogenes',
	'Clostridium_thermosulfurogenes', 'Desulfotomaculum_alkaliphilum',
	'Enterococcus_faecalis', 'Geobacillus_thermodenitrificans',
	'Geobacillus_thermoleovorans', 'Halonatronum_saccharophilum',
	'Lactobacillus_acidophilus', 'Lactobacillus_delbrueckii',
	'Lactobacillus_paracasei', 'Lactobacillus_rhamnosus',
	'Lactococcus_lactis', 'Listeria_monocytogenes', 'Pelotomaculum_thermopropionicum',
	'Planococcus_halocryophilus', 'Staphylococcus_aureus', 'Staphylococcus_xylosus',
	'Streptococcus_salivarius', 'Streptococcus_thermophilus', 'Sulfobacillus_sibiricus',
	'Sulfobacillus_thermotolerans', 
	'Thermacetogenium_phaeum', 'Thermaerobacter_nagasakiensis',
	'Thermoanaerobacter_ethanolicus', 'Thermoanaerobacter_keratinophilus',
	'Thermoanaerobacter_mathranii', 'Thermoanaerobacter_siderophilus',
	'Thermoanaerobacter_subterraneus', 'Thermoanaerobacter_tengcongensis',
	'Thermoanaerobacter_yonseiensis', 'Thermobrachium_celere', 
	'Trichococcus_patagoniensis'
)

Proteobacteria <- c(
	'Acidithiobacillus_ferrivorans', 'Acidocella_aromatica',
	'Aeromonas_hydrophila', 'Colwellia_demingiae', 'Colwellia_hornerae',
	'Colwellia_psychrerythraea', 'Colwellia_psychrotropica', 
	'Desulfobacter_curvatus', 'Desulfofaba_gelida', 'Desulfofrigus_fragile',
	'Desulfofrigus_oceanense', 'Desulforhopalus_vacuolatus',
	'Desulfotalea_arctica', 'Desulfotalea_psychrophila',
	'Erwinia_amylovora', 'Escherichia_coli', 'Glaciecola_punicea',
	'Halomonas_campisalis', 'Halomonas_elongata', 'Halomonas_subglaciescola',
	'Hydrogenophaga_pseudoflava', 'Hydrogenophilus_hirschii',
	'Klebsiella_oxytoca', 'Marinobacter_alkaliphilus',
	'Paracoccus_halodenitrificans', 'Pseudoalteromonas_antarctica',
	'Pseudoalteromonas_haloplanktis', 'Pseudomonas_aeruginosa', 
	'Pseudomonas_fluorescens', 'Pseudomonas_putida', 'Psychrobacter_glacincola',
	'Psychrobacter_muriicola', 'Serratia_marcescens', 'Shewanella_gelidimarina',
	'Vibrio_marinus', 'Yersinia_enterocolitica'
)

# Find the species from these 5 phyla and remove the other species.
all_dat_E$Phylum <- NA
all_dat_E$Phylum[all_dat_E$Species %in% Bacillariophyta] <- 'Bacillariophyta'
all_dat_E$Phylum[all_dat_E$Species %in% Cyanobacteria] <- 'Cyanobacteria'
all_dat_E$Phylum[all_dat_E$Species %in% Dinophyta] <- 'Dinophyta'
all_dat_E$Phylum[all_dat_E$Species %in% Proteobacteria] <- 'Proteobacteria'
all_dat_E$Phylum[all_dat_E$Species %in% Firmicutes] <- 'Firmicutes'
all_dat_E$Phylum[all_dat_E$Species %in% Euryarchaeota] <- 'Euryarchaeota'
all_dat_E <- all_dat_E[!is.na(all_dat_E$Phylum),]

all_dat_W_op$Phylum <- NA
all_dat_W_op$Phylum[all_dat_W_op$Species %in% Bacillariophyta] <- 'Bacillariophyta'
all_dat_W_op$Phylum[all_dat_W_op$Species %in% Cyanobacteria] <- 'Cyanobacteria'
all_dat_W_op$Phylum[all_dat_W_op$Species %in% Dinophyta] <- 'Dinophyta'
all_dat_W_op$Phylum[all_dat_W_op$Species %in% Proteobacteria] <- 'Proteobacteria'
all_dat_W_op$Phylum[all_dat_W_op$Species %in% Firmicutes] <- 'Firmicutes'
all_dat_W_op$Phylum[all_dat_W_op$Species %in% Euryarchaeota] <- 'Euryarchaeota'
all_dat_W_op <- all_dat_W_op[!is.na(all_dat_W_op$Phylum),]

# Plot the resulting distributions.
p1 <- ggplot(all_dat_E, aes(x = exp(ln_E), color = Phylum)) +
	geom_density(alpha = 0.3, lwd = 1.2) +
	xlim(c(0,4)) +
	scale_color_manual(
		values = c(
			'Proteobacteria' = '#af9fcc',
			'Bacillariophyta' = '#bd8891', 
			'Cyanobacteria' = '#7fd7e9', 
			'Dinophyta' = '#7fc0d9', 
			'Firmicutes' = '#d1e9c7', 
			'Euryarchaeota' = '#ffc7c7'
		)
	) +
	guides(fill = FALSE) +
	xlab(expression(italic(E) * " (eV)")) +
	theme(
		axis.text = element_text(size = 8),
		axis.title = element_text(size = 10),
		legend.title = element_text(size = 10, face = 'bold'),
		legend.text = element_text(size = 9),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center'
	)
	
ggsave(p1, file = '../Results/E_distributions.pdf', width = 4.5, height = 3, units = 'in')

p2 <- ggplot(all_dat_W_op, aes(x = exp(ln_W_op), color = Phylum)) +
	geom_density(alpha = 0.3, lwd = 1.2) +
	xlim(c(3,32.5)) +
	scale_color_manual(
		values = c(
			'Proteobacteria' = '#af9fcc',
			'Bacillariophyta' = '#bd8891', 
			'Cyanobacteria' = '#7fd7e9', 
			'Dinophyta' = '#7fc0d9', 
			'Firmicutes' = '#d1e9c7', 
			'Euryarchaeota' = '#ffc7c7'
		)
	) +
	guides(fill = FALSE) +
	xlab(expression(italic(W)[op] * " (°C)")) +
	theme(
		axis.text = element_text(size = 8),
		axis.title = element_text(size = 10),
		legend.title = element_text(size = 10, face = 'bold'),
		legend.text = element_text(size = 9),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center'
	)
	
ggsave(p2, file = '../Results/W_op_distributions.pdf', width = 4.5, height = 3, units = 'in')

p_Proteobacteria_E <- ggplot(all_dat_E[all_dat_E$Phylum == "Proteobacteria",], aes(x = exp(ln_E))) + 
	ylim(c(0, 1.15)) +
	xlim(c(0, 4)) +
	xlab(expression(italic(E) * " (eV)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#af9fcc', fill = '#af9fcc') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 6),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 3.2, y = 1, label = 'italic(N) == 36', parse = TRUE)

p_Proteobacteria_W_op <- ggplot(all_dat_W_op[all_dat_W_op$Phylum == "Proteobacteria",], aes(x = exp(ln_W_op))) + 
	ylim(c(0, 0.15)) +
	xlim(c(3, 32.5)) +
	ylab('') +
	xlab(expression(italic(W)[op] * " (°C)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#af9fcc', fill = '#af9fcc') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 9),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 28, y = 0.12, label = 'italic(N) == 32', parse = TRUE)

ggsave(
	plot_grid(p_Proteobacteria_E, p_Proteobacteria_W_op, ncol = 2), 
	file = '../Results/Proteobacteria.pdf',
	width = 5, height = 2, units = 'in'
)

p_Firmicutes_E <- ggplot(all_dat_E[all_dat_E$Phylum == "Firmicutes",], aes(x = exp(ln_E))) + 
	ylim(c(0, 1.15)) +
	xlim(c(0, 4)) +
	xlab(expression(italic(E) * " (eV)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#d1e9c7', fill = '#d1e9c7') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 6),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 3.2, y = 1, label = 'italic(N) == 52', parse = TRUE)

p_Firmicutes_W_op <- ggplot(all_dat_W_op[all_dat_W_op$Phylum == "Firmicutes",], aes(x = exp(ln_W_op))) + 
	ylim(c(0, 0.15)) +
	xlim(c(3, 32.5)) +
	ylab('') +
	xlab(expression(italic(W)[op] * " (°C)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#d1e9c7', fill = '#d1e9c7') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 9),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 28, y = 0.12, label = 'italic(N) == 40', parse = TRUE)

ggsave(
	plot_grid(p_Firmicutes_E, p_Firmicutes_W_op, ncol = 2), 
	file = '../Results/Firmicutes.pdf',
	width = 5, height = 2, units = 'in'
)

p_Cyanobacteria_E <- ggplot(all_dat_E[all_dat_E$Phylum == "Cyanobacteria",], aes(x = exp(ln_E))) + 
	ylim(c(0, 1.15)) +
	xlim(c(0, 4)) +
	xlab(expression(italic(E) * " (eV)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#7fd7e9', fill = '#7fd7e9') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 6),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 3.2, y = 1, label = 'italic(N) == 17', parse = TRUE)

p_Cyanobacteria_W_op <- ggplot(all_dat_W_op[all_dat_W_op$Phylum == "Cyanobacteria",], aes(x = exp(ln_W_op))) + 
	ylim(c(0, 0.15)) +
	xlim(c(3, 32.5)) +
	ylab('') +
	xlab(expression(italic(W)[op] * " (°C)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#7fd7e9', fill = '#7fd7e9') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 9),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 28, y = 0.12, label = 'italic(N) == 11', parse = TRUE)

ggsave(
	plot_grid(p_Cyanobacteria_E, p_Cyanobacteria_W_op, ncol = 2), 
	file = '../Results/Cyanobacteria.pdf',
	width = 5, height = 2, units = 'in'
)

p_Euryarchaeota_E <- ggplot(all_dat_E[all_dat_E$Phylum == "Euryarchaeota",], aes(x = exp(ln_E))) + 
	ylim(c(0, 1.15)) +
	xlim(c(0, 4)) +
	xlab(expression(italic(E) * " (eV)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#ffc7c7', fill = '#ffc7c7') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 6),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 3.2, y = 1, label = 'italic(N) == 33', parse = TRUE)

p_Euryarchaeota_W_op <- ggplot(all_dat_W_op[all_dat_W_op$Phylum == "Euryarchaeota",], aes(x = exp(ln_W_op))) + 
	ylim(c(0, 0.15)) +
	xlim(c(3, 32.5)) +
	ylab('') +
	xlab(expression(italic(W)[op] * " (°C)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#ffc7c7', fill = '#ffc7c7') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 9),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 28, y = 0.12, label = 'italic(N) == 27', parse = TRUE)

ggsave(
	plot_grid(p_Euryarchaeota_E, p_Euryarchaeota_W_op, ncol = 2), 
	file = '../Results/Euryarchaeota.pdf',
	width = 5, height = 2, units = 'in'
)

p_Bacillariophyta_E <- ggplot(all_dat_E[all_dat_E$Phylum == "Bacillariophyta",], aes(x = exp(ln_E))) + 
	ylim(c(0, 1.15)) +
	xlim(c(0, 4)) +
	xlab(expression(italic(E) * " (eV)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#bd8891', fill = '#bd8891') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 6),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 3.2, y = 1, label = 'italic(N) == 27', parse = TRUE)

p_Bacillariophyta_W_op <- ggplot(all_dat_W_op[all_dat_W_op$Phylum == "Bacillariophyta",], aes(x = exp(ln_W_op))) + 
	ylim(c(0, 0.15)) +
	xlim(c(3, 32.5)) +
	ylab('') +
	xlab(expression(italic(W)[op] * " (°C)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#bd8891', fill = '#bd8891') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 9),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 28, y = 0.12, label = 'italic(N) == 14', parse = TRUE)

ggsave(
	plot_grid(p_Bacillariophyta_E, p_Bacillariophyta_W_op, ncol = 2), 
	file = '../Results/Bacillariophyta.pdf',
	width = 5, height = 2, units = 'in'
)

p_Dinophyta_E <- ggplot(all_dat_E[all_dat_E$Phylum == "Dinophyta",], aes(x = exp(ln_E))) + 
	ylim(c(0, 1.15)) +
	xlim(c(0, 4)) +
	xlab(expression(italic(E) * " (eV)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#7fc0d9', fill = '#7fc0d9') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 6),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 3.2, y = 1, label = 'italic(N) == 19', parse = TRUE)

p_Dinophyta_W_op <- ggplot(all_dat_W_op[all_dat_W_op$Phylum == "Dinophyta",], aes(x = exp(ln_W_op))) + 
	ylim(c(0, 0.15)) +
	xlim(c(3, 32.5)) +
	ylab('') +
	xlab(expression(italic(W)[op] * " (°C)")) +
	geom_density(alpha = 0.3, lwd = 1.2, color = '#7fc0d9', fill = '#7fc0d9') +
	theme(
		axis.text = element_text(size = 6*1.5),
		axis.title = element_text(size = 8*1.5),
		legend.title = element_text(size = 6, face = 'bold'),
		legend.text = element_text(size = 9),
		plot.margin=unit(c(1,1,0,0),"mm"),
		legend.margin=margin(0,0,0,0),
		legend.box.margin=margin(-5,0,0,0),
		legend.justification = 'center',
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
	) +
	annotate('text', x = 28, y = 0.12, label = 'italic(N) == 15', parse = TRUE)

ggsave(
	plot_grid(p_Dinophyta_E, p_Dinophyta_W_op, ncol = 2), 
	file = '../Results/Dinophyta.pdf',
	width = 5, height = 2, units = 'in'
)
