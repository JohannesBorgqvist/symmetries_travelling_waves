#=================================================================================
#=================================================================================
# Script:"plot_travelling_wave_symmetries"
# Date: 2022-07-22
# Implemented by: Johannes Borgqvist
# Description:
# The script plots all the travelling wave symmetries of the generalised growth model and the growth model with an Allee effect.
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
from numpy import *
import matplotlib.pyplot as plt
#=================================================================================
#=================================================================================
# Functions
#=================================================================================
#=================================================================================
# This functions enable us to reproduce our plots using pgfplots in LaTeX
def plot_LaTeX_2D(t,y,file_str,plot_str,legend_str):
    # Open a file with the append option
    # so that we can write to the same
    # file multiple times
    f = open(file_str, "a")
    # Create a temporary string which
    # is the one that does the plotting.
    # Here we incorporate the input plot_str
    # which contains the color, and the markers
    # of the plot at hand
    if len(legend_str)==0:
        temp_str = "\\addplot[\nforget plot,\n" + plot_str+ "\n]\n"
    else:
        temp_str = "\\addplot[\n" + plot_str+ "\n]\n"
    # Add the coordinates
    temp_str += "coordinates {%\n"
    # Loop over the input files and add
    # them to the file
    for i in range(len(t)):
        temp_str += "(" + str(t[i]) + "," + str(y[i]) + ")\n"
    # The plotting is done, let's close the shop    
    temp_str += "};\n"
    # Add a legend if one is provided
    if len(legend_str) > 0:
        temp_str += "\\addlegendentry{" + legend_str + "}\n"
    # Finally, we write the huge string
    # we have created
    f.write("%s"%(temp_str))
    # Close the file
    f.close()
#=================================================================================
#=================================================================================
# Plot the travelling wave solutions
#=================================================================================
#=================================================================================
# Define the overall parameters
C= 0.5
A = 0.5
B = 0.5
kappa = 5
# Define our z vector
z = linspace(0,20,500)

# Define our travelling wave for the Allee model
# Allee model
lambda_allee = sqrt(kappa/2)
tv_allee = array([((C*exp(lambda_allee*z_temp))/(1+C*exp(lambda_allee*z_temp))) for z_temp in z])
# Plot these solutions
f1, ax_1 = plt.subplots(1, 1, constrained_layout=True, figsize=(20, 8))
ax_1.plot(z, tv_allee, '-', label="Allee model" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_1.grid()
ax_1.legend(loc='best',prop={"size":20})
# Set fontsize of labels
ax_1.set_xlabel(xlabel='Wave front, $z=x-ct$',fontsize=25)
ax_1.set_ylabel(ylabel='Travelling wave, $u(z)$',fontsize=25)
# Change the size of the ticks
ax_1.tick_params(axis='both', which='major', labelsize=20)
ax_1.tick_params(axis='both', which='minor', labelsize=20)
# Title and saving the figure
f1.suptitle('Travelling wave solution of RD-models of cell growth',fontsize=30,weight='bold')
f1.savefig('../Figures/travelling_wave_solution_Allee.png')
# Show the plot at last
plt.show()


# STABLE POINT RADIAL SYMMETRY
plot_LaTeX_2D(z, tv_allee,"../Figures/travelling_waves/Input/tv_1.tex","color=clr_1,line width=1.5pt,","Allee model")

