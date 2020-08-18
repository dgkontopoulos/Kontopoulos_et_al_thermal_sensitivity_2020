#!/usr/bin/Rscript

# This script performs a numerical exploration of the relationship 
# between E and W_op in the Sharpe-Schoolfield model (S1 Appendix, Fig 
# Ea).
#
# It also generates a plot of three TPCs with common values for B_0, 
# T_pk, and E_D, but different values for E (S1 Appendix, Fig Eb).

library(cowplot)
library(R.devices)

#####################
# F U N C T I O N S #
#####################

# This function calculates W_op, given trait and temperature values.
find_W_op <- function(TPC, temps)
{
	
	# Find the temperature at which the TPC peaks (T_pk), and the value 
	# at the peak (B_pk).
	T_pk <- temps[TPC == max(TPC)]
	B_pk <- max(TPC)
	
	
	# Find the temperature at the rise of the TPC where trait performance 
	# is ~ half of B_pk.
	min_difference <- 9999
	current_temp <- NA
	
	for ( i in 1:length(temps) )
	{
		if ( temps[i] >= T_pk )
		{
			break
		}
		
		if ( sqrt((TPC[i] - (B_pk/2))^2) < min_difference )
		{
			current_temp <- temps[i]
			min_difference <- sqrt((TPC[i] - (B_pk/2))^2)
		}
	}
	
	# Calculate W_op as the difference between T_pk and the temperature 
	# at the rise of the TPC where trait performance is ~ 0.5 times B_pk.
	return(T_pk - current_temp)
	
}

# This function calculates trait performance based on a 4-parameter 
# variant of the Sharpe-Schoolfield model.
#
# T_ref is set to 273.15 K.
schoolf <- function(B0, E, E_D, T_pk, temp)
{
	
	# Set the value of the Boltzmann constant (k).
	k <- 8.617 * 10^-5
	
	return(B0 * exp(-E * ((1/(k*temp)) - (1/(k*273.15))))/(1 + (E/(E_D - E)) * exp(E_D/k * (1/T_pk - 1/temp))))
}

####################
# M A I N  C O D E #
####################

################
# P a n e l  A #
################

# Get a vector of temperatures from 0 to 60 degrees Celsius.
temps <- seq(0,60,0.05) + 273.15

# Get a vector of E values from 0.15 to 3.95 eV.
Es <- seq(0.15, 3.95, 0.05)

# Prepare lists and vectors to store trait values and W_op values.
TPC1s <- list()
TPC2s <- list()
TPC3s <- list()

W_op1s <- c()
W_op2s <- c()
W_op3s <- c()

counter <- 1

# For each E value...
for ( i in Es )
{

	# ... calculate trait values for 3 parameter combinations...
	TPC1s[[counter]] <- schoolf(5e-3, i, 5.5, 30 + 273.15, temps)
	TPC2s[[counter]] <- schoolf(1e-3, i, 5.5, 40 + 273.15, temps)
	TPC3s[[counter]] <- schoolf(1e-3, i, 7.5, 30 + 273.15, temps)
	
	# ... and calculate W_op values.
	W_op1s[counter] <- find_W_op(TPC1s[[counter]], temps)
	W_op2s[counter] <- find_W_op(TPC2s[[counter]], temps)
	W_op3s[counter] <- find_W_op(TPC3s[[counter]], temps)
	
	counter <- counter + 1
}

# Plot the resulting relationships between E and W_op.
dat_to_plot <- data.frame(
	E = rep(Es, 3),
	W_op = c(W_op1s, W_op2s, W_op3s),
	TPC = c(
		rep('1', length(Es)), rep('2', length(Es)), 
		rep('3', length(Es))
	)
)

p1 <- ggplot(dat_to_plot, aes(x = W_op, y = E, color = as.factor(TPC))) + 
	geom_line(lwd = 1) +
	ylab('Rate value') +
	ylim(c(0,4)) +
	ylab(bquote(italic(E) ~ '(eV)')) +
	xlab(bquote(italic(W)[op] ~ '(°C)')) +
	xlim(c(min(dat_to_plot$W_op), 38)) +
	scale_color_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb')) +
	guides(color = FALSE) +
	theme(
		axis.text.y = element_text(size = 8),
		axis.text.x = element_text(size = 8),
		plot.title = element_text(size = 10),
		axis.title = element_text(size = 10),
		legend.title = element_text(size = 10, face = 'bold'),
		legend.text = element_text(size = 10),
		plot.margin=unit(c(0,0,0,0),"mm")
	)

suppressGraphics(
	ggsave(
		p1, file = '../Results/E_vs_W_op_with_other_parameters_changing.pdf', 
		width = 4.5 * 1.5 / 2, height = 2 * 1.25, units = 'in'
	)
)

################
# P a n e l  B #
################

# Get a vector of temperatures from -5 to 60 degrees Celsius.
temps <- seq(-5,60,0.25) + 273.15

# Get three TPCs that have the same B_0, T_pk, and E_D values, but 
# different E values.
TPC_1 <- schoolf(1e-3, 0.3, 1.5, 30 + 273.15, temps)
TPC_2 <- schoolf(1e-3, 0.6, 1.5, 30 + 273.15, temps)
TPC_3 <- schoolf(1e-3, 1, 1.5, 30 + 273.15, temps)

# Print the W_op values for the 3 TPCs.
print(find_W_op(TPC_1, temps))
print(find_W_op(TPC_2, temps))
print(find_W_op(TPC_3, temps))

# Prepare a data frame for plotting the TPCs.
dat_to_plot <- data.frame(
	temps = rep(temps, 3),
	values = c(TPC_1, TPC_2, TPC_3),
	fit = c(
		rep('fit_1', length(temps)), rep('fit_2', length(temps)),
		rep('fit_3', length(temps))
	)
)

# Prepare lines that start from T_pk and extend to the temperature where 
# trait performance is half of B_pk.
line_W_op_1 <- data.frame(
	x = seq((30 + 273.15) - find_W_op(TPC_1, temps), 30 + 273.15, 0.05), 
	y = rep(max(TPC_1)/2, length(seq((30 + 273.15) - find_W_op(TPC_1, temps), 30 + 273.15, 0.05)))
)

line_W_op_2 <- data.frame(
	x = seq((30 + 273.15) - find_W_op(TPC_2, temps), 30 + 273.15, 0.05), 
	y = rep(max(TPC_2)/2, length(seq((30 + 273.15) - find_W_op(TPC_2, temps), 30 + 273.15, 0.05)))
)

line_W_op_3 <- data.frame(
	x = seq((30 + 273.15) - find_W_op(TPC_3, temps), 30 + 273.15, 0.05), 
	y = rep(max(TPC_3)/2, length(seq((30 + 273.15) - find_W_op(TPC_3, temps), 30 + 273.15, 0.05)))
)

# Plot the three TPCs.
p2 <- ggplot(dat_to_plot, aes(x = temps - 273.15, y = values, color = fit)) + 
	geom_line(lwd = 1.3) +
	geom_line(data = line_W_op_1, aes(x = x - 273.15, y = y), lwd = 1, lty = 1, color = '#66c2a5') +
	geom_line(data = line_W_op_2, aes(x = x - 273.15, y = y), lwd = 1, lty = 1, color = '#fc8d62') +
	geom_line(data = line_W_op_3, aes(x = x - 273.15, y = y), lwd = 1, lty = 1, color = '#8da0cb') +
	ylab('Rate value') +
	xlab('Temperature (°C)') +
	xlim(c(-5,50)) +
	geom_vline(xintercept = 30, lty = 2, lwd = 0.5) + 
	scale_color_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb')) +
	guides(color = FALSE) +
	ylim(c(0, max(dat_to_plot$values))) +
	theme(
		axis.text.y = element_text(size = 8),
		axis.text.x = element_text(size = 8),
		plot.title = element_text(size = 10),
		axis.title = element_text(size = 10),
		legend.title = element_text(size = 10, face = 'bold'),
		legend.text = element_text(size = 10),
		plot.margin=unit(c(0,0,0,0),"mm")
	)

suppressGraphics(
	ggsave(
		p2, file = '../Results/TPCs_with_only_different_E.pdf', 
		width = 4.5 * 1.5 / 2, height = 2 * 1.25, units = 'in'
	)
)
