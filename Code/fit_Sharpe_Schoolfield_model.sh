#!/bin/sh

# This script calls the TPC_fitting.py script to fit the Sharpe-Schoolfield 
# model to the 4 datasets of this study and obtain estimates of TPC 
# parameters.

/usr/bin/python TPC_fitting.py ../Data/phytoplankton_specific_growth_rate.csv ../Data/TPC_parameter_estimates_phytoplankton_r_max.csv
/usr/bin/python TPC_fitting.py ../Data/prokaryotes_specific_growth_rate.csv ../Data/TPC_parameter_estimates_prokaryotes_r_max.csv
/usr/bin/python TPC_fitting.py ../Data/photosynthesis_rate.csv ../Data/TPC_parameter_estimates_plants_net_photosynthesis_rate.csv
/usr/bin/python TPC_fitting.py ../Data/respiration_rate.csv ../Data/TPC_parameter_estimates_plants_respiration_rate.csv
