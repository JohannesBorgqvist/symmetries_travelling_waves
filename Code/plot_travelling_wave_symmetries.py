#=================================================================================
#=================================================================================
# Script:"plot_travelling_wave_symmetries"
# Date: 2023-02-02
# Implemented by: Johannes Borgqvist
# Description:
# The script plots all the travelling wave symmetries of the Fisher--KPP model
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
from numpy import * # For arrays and numerics
import matplotlib.pyplot as plt # For plotting
# Scipy to do all the computations and integrations
from scipy import integrate # For solving ODEs.
# Make the plotting fancy with LaTeX
plt.rcParams['text.usetex'] = True
#=================================================================================
#=================================================================================
# Functions
#=================================================================================
#=================================================================================
# Function 1: dX_dt_Fisher_KPP
# The travelling wave Fisher KPP model in terms of
# of the pseudo state v=u'. 
def dX_dt_Fisher_KPP(X, t=0,*parameters):
    # Return the dynamics of the logistic model
    return array([ X[1],
                   -(5/sqrt(6))*X[1]-X[0]*(1-X[0])])
# Function 2: Gamma_Feng
# This function plots the action of Feng's symmetry.
# The first component of this ODE system is xi, the second
# is eta and the third one is eta^{(1)}.
def Gamma_Feng(X, epsilon=0):
    # Define c
    c = 5/sqrt(6)
    # Return the dynamics of the logistic symmetry
    return array([ -exp((c*X[0])/(5))*((5)/(2*c)),
                  exp((c*X[0])/(5))*X[1],
                exp((c*X[0])/(5))*((c/(5))*X[1]+(3/2)*X[2])])
# Function 3: For plotting in LaTeX
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
# Solving the travelling wave Fisher KPP
#=================================================================================
#=================================================================================
# Define the initial conditions of the original solution
X_0_Feng = array([0.75, 0])
# Define the initial conditions of the original solution
X_0_trans = array([0.95, 0])
# Define the travelling wave vector
z_vec = linspace(0,15,300)
# Solve the ODE for the travelling wave equation
X_Fisher, infodict = integrate.odeint(dX_dt_Fisher_KPP, X_0_trans, z_vec,full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the solution to our lovely system
u_Fisher_1, v_Fisher_1 = X_Fisher.T
# Solve the ODE for the travelling wave equation
X_Fisher, infodict = integrate.odeint(dX_dt_Fisher_KPP, X_0_Feng, z_vec,full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the solution to our lovely system
u_Fisher_2, v_Fisher_2 = X_Fisher.T
#=================================================================================
#=================================================================================
# Transforming solutions of the travelling wave Fisher-KPP model with the
# translation symmetry
#=================================================================================
#=================================================================================
# Define some magical indices
magical_indices = arange(0,len(z_vec),20)
# Epsilon
epsilon = 2.0
epsilon_vec = linspace(0,epsilon,50)
# Allocate memory for the transformed states
z_trans = []
u_trans = []
v_trans = []
# Loop over the magical indices
for magical_index in magical_indices:
    # Solve the ODE for the symmetry
    z_temp = array([z_vec[magical_index]+epsilon_temp for epsilon_temp in epsilon_vec])
    u_temp = array([u_Fisher_1[magical_index]+epsilon_temp-epsilon_temp for epsilon_temp in epsilon_vec])
    v_temp = array([v_Fisher_1[magical_index]+epsilon_temp-epsilon_temp for epsilon_temp in epsilon_vec])    
    # Save our solutions here
    z_trans.append(z_temp)
    u_trans.append(u_temp)
    v_trans.append(v_temp)
# Define a transformed z vector
z_vec_trans = linspace(z_trans[0][-1],z_vec[-1],300)
# Plot the transformed solution
X_Fisher_trans, infodict = integrate.odeint(dX_dt_Fisher_KPP, array([u_trans[0][-1], v_trans[0][-1]]), z_vec_trans,full_output=True)
# Extract the solution to our lovely system
u_Fisher_trans, v_Fisher_trans = X_Fisher_trans.T
#=================================================================================
#=================================================================================
# Transforming solutions of the travelling wave Fisher-KPP model with Feng's
# symmetry
#=================================================================================
#=================================================================================
# Define some magical indices
magical_indices = arange(0,round(4*len(z_vec)/5),20)
# Epsilon
epsilon = 0.25
epsilon_vec = linspace(0,epsilon,50)
# Allocate memory for the transformed states
z_Feng = []
u_Feng = []
v_Feng = []
# Loop over the magical indices
for magical_index in magical_indices:
    # Calculate the action of Feng's symmetry by solving the ODE system
    X_Fisher_Feng, infodict = integrate.odeint(Gamma_Feng, array([z_vec[magical_index],u_Fisher_2[magical_index], v_Fisher_2[magical_index]]), epsilon_vec,full_output=True)
    # Extract the solution to our lovely system
    z_temp,u_temp, v_temp = X_Fisher_Feng.T
    # Save our solutions here
    z_Feng.append(z_temp)
    u_Feng.append(u_temp)
    v_Feng.append(v_temp)
# Define a transformed z vector
z_vec_Feng = linspace(z_Feng[0][-1],z_vec[-1],300)
# Plot the transformed solution
X_Fisher_Feng, infodict = integrate.odeint(dX_dt_Fisher_KPP, array([u_Feng[0][-1], v_Feng[0][-1]]), z_vec_Feng,full_output=True)
# Extract the solution to our lovely system
u_Fisher_Feng, v_Fisher_Feng = X_Fisher_Feng.T
#=================================================================================
#=================================================================================
# Plotting the symmetries of the logistic model
#=================================================================================
#=================================================================================
#Define the first figure
f_sym_tv, ax_sym_tv = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
#---------------------------------------------------------------------------------
# Subplot 1 of 2: The translation symmetry
# The original solution
ax_sym_tv[0].plot(z_vec,u_Fisher_1,color=(0/256,68/256,27/256),label="$u(z)$",linewidth=1.0)
# The symmetry
for index in range(len(magical_indices)):
    if index == 0:
        ax_sym_tv[0].plot(z_trans[index], u_trans[index], '--', label="$\\Gamma_{\\epsilon}^{z}$",color=(0/256,0/256,0/256),linewidth=1.0)
    else:
        ax_sym_tv[0].plot(z_trans[index], u_trans[index], '--',color=(0/256,0/256,0/256),linewidth=1.0)
# The transformed solution        
ax_sym_tv[0].plot(z_vec_trans,u_Fisher_trans,color=(153/256,216/256,201/256),label="$\\hat{u}(z)$",linewidth=1.0)
# Set the y-limit
ax_sym_tv[0].set_ylim([0, 1])
# Set a grid and define a legend
ax_sym_tv[0].grid()
ax_sym_tv[0].legend(loc='best',prop={"size":15})
# Set the x-labels and y-labels
ax_sym_tv[0].set_xlabel(xlabel="Travelling wave variable, $z=x-ct$",fontsize=15)
ax_sym_tv[0].set_ylabel(ylabel="Population density, $u(z)$",fontsize=15)
ax_sym_tv[0].set_title(label="Translation symmetry",fontsize=20)
#---------------------------------------------------------------------------------
# Subplot 1 of 2: Feng's symmetry
# The original solution
ax_sym_tv[1].plot(z_vec,u_Fisher_2,color=(0/256,68/256,27/256),label="$u(z)$",linewidth=1.0)
# The symmetry
for index in range(len(magical_indices)):
    if index == 0:
        ax_sym_tv[1].plot(z_Feng[index], u_Feng[index], '--', label="$\\Gamma_{\\epsilon}^{F}$",color=(0/256,0/256,0/256),linewidth=1.0)
    else:
        ax_sym_tv[1].plot(z_Feng[index], u_Feng[index], '--',color=(0/256,0/256,0/256),linewidth=1.0)
# The transformed solution        
ax_sym_tv[1].plot(z_vec_Feng,u_Fisher_Feng,color=(153/256,216/256,201/256),label="$\\hat{u}(z)$",linewidth=1.0)
# Set the y-limit
ax_sym_tv[1].set_ylim([0, 1])
# Set a grid and define a legend
ax_sym_tv[1].grid()
ax_sym_tv[1].legend(loc='best',prop={"size":15})
# Set the x-labels and y-labels
ax_sym_tv[1].set_xlabel(xlabel="Travelling wave variable, $z=x-ct$",fontsize=15)
ax_sym_tv[1].set_ylabel(ylabel="Population density, $u(z)$",fontsize=15)
ax_sym_tv[1].set_title(label="Feng's symmetry",fontsize=20)
plt.show()
#=================================================================================
#=================================================================================
# SAVE SOLUTIONS IN LATEX AS WELL
#=================================================================================
#=================================================================================
# THE TRANSLATION SYMMETRY
plot_LaTeX_2D(z_vec,u_Fisher_1,"../Figures/travelling_waves/Input/translation.tex","color=log_1,line width=1pt,","$u(z)$")
plot_LaTeX_2D(z_vec_trans,u_Fisher_trans,"../Figures/travelling_waves/Input/translation.tex","color=log_2,line width=1pt,","$\\hat{u}(z)$")
# The symmetry
for index in range(len(magical_indices)):
    if index == 0:
        plot_LaTeX_2D(z_trans[index], u_trans[index],"../Figures/travelling_waves/Input/translation.tex","color=black,->,>=latex,densely dashed","$\\Gamma_{\\epsilon}^{Z}$")
    else:
        plot_LaTeX_2D(z_trans[index], u_trans[index],"../Figures/travelling_waves/Input/translation.tex","color=black,->,>=latex,densely dashed",[])
# THE FENG'S SYMMETRY
plot_LaTeX_2D(z_vec,u_Fisher_2,"../Figures/travelling_waves/Input/Feng.tex","color=log_1,line width=1pt,","$u(z)$")
plot_LaTeX_2D(z_vec_Feng,u_Fisher_Feng,"../Figures/travelling_waves/Input/Feng.tex","color=log_2,line width=1pt,","$\\hat{u}(z)$")
# The symmetry
for index in range(len(magical_indices)):
    if index == 0:
        plot_LaTeX_2D(z_Feng[index], u_Feng[index],"../Figures/travelling_waves/Input/Feng.tex","color=black,->,>=latex,densely dashed","$\\Gamma_{\\epsilon}^{F}$")
    else:
        plot_LaTeX_2D(z_Feng[index], u_Feng[index],"../Figures/travelling_waves/Input/Feng.tex","color=black,->,>=latex,densely dashed",[])
