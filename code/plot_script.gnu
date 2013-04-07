#
# This script takes data files from "ising.f90" and plots them
#

# Set terminal to create png output

set term png



#######################################################
### Magnetization vs. Time for Metropolis Algorithm ###
#######################################################

set ylabel 'Magnetization (arbitrary units)'
set xlabel 'Time (arbitrary units)'

set output 'metrop_mag_time.png'
set title 'Magnetization vs. Time -- (Metropolis)'
plot 'metrop_mag_time.txt' using 2:1



#######################################
#### Magnetization vs. Temperature ####
#######################################


set xlabel 'Temperature (arbitrary units)'

set output 'metrop_mag_temp.png'
set title 'Magnetization vs. Temperature -- (Metropolis)'
plot 'metrop_mag_temp.txt' using 2:1


set output 'wolff_mag_temp.png'
set title 'Magnetization vs. Temperature -- (Wolff)'
plot 'wolff_mag_temp.txt' using 2:1


set output 'swenwang_mag_temp.png'
set title 'Magnetization vs. Temperature -- (Swensden-Wang)'
plot 'swenwang_mag_temp.txt' using 2:1



########################################
### Magnetic Susceptibility vs. Temp ###
########################################

set ylabel 'Magnetic Susceptibility (arbitrary units)'

set output 'metrop_magsus_temp.png'
set title 'Magnetic Susceptibility vs. Temperature -- (Metropolis)'
plot 'metrop_magsus_temp.txt' using 2:1

set output 'wolff_magsus_temp.png'
set title 'Magnetic Susceptibility vs. Temperature -- (Wolff)'
plot 'wolff_magsus_temp.txt' using 2:1


set output 'swenwang_magsus_temp.png'
set title 'Magnetic Susceptibility vs. Temperature -- (Swensden-Wang)'
plot 'swenwang_magsus_temp.txt' using 2:1


######################
####### Energy #######
######################

#set yrange [3e+10: 5e+10]

#set output 'metrop_energy.png'
#set title 'Energy vs. Temperature -- (Metropolis)'
#plot 'metrop_energy.txt' using 2:1

#set output 'wolff_energy.png'
#set title 'Energy vs. Temperature -- (Wolff)'
#plot 'wolff_energy.txt' using 2:1

#set output 'swen_energy.png'
#set title 'Energy vs. Temperature -- (Swensden-Wang)'
#plot 'swen_energy.txt' using 2:1


