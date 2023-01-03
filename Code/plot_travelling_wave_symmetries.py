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
# Scipy to do all the computations and integrations
from scipy import integrate
from scipy.optimize import fsolve
from scipy import integrate # For solving ODEs.
# Make the plotting fancy with LaTeX
plt.rcParams['text.usetex'] = True
#=================================================================================
#=================================================================================
# Functions
#=================================================================================
#=================================================================================
# Function 1: dX_dt_exponential
# The second order exponential model as a system of two ODEs in terms
# of the pseudo state v=u'. 
def dX_dt_exponential(X, t=0,*parameters):
    # Extract the wave speed c
    c = parameters[0]    
    # Return the dynamics of the exponential model
    return array([ X[1],
                   -c*X[1]-X[0]])
# Function 2: Gamma_exponential
# This function plots the symmetry of the exponential model. The first
# Component of this ODE system is xi=1, the second is eta=u and the
# third one is eta^{(1)}=u'.
def Gamma_exponential(X, epsilon=0):  
    # Return the dynamics of the exponential symmetry
    return array([ 1,
                   X[1],
                   X[2]])
# Function 3: dX_dt_logistic
# The second order logistic model as a system of two ODEs in terms
# of the pseudo state v=u'. 
def dX_dt_logistic(X, t=0,*parameters):
    # Extract the travelling wave parameter
    c = parameters[0]
    k1 = parameters[1]
    # Return the dynamics of the logistic model
    return array([ X[1],
                   -(c-k1)*X[1]-X[0]*(1-X[0])])
# Function 4: Gamma_logistic
# This function plots the symmetry of the logistic model. The first
# Component of this ODE system is xi, the second is eta and the
# third one is eta^{(1)}.
def Gamma_logistic(X, epsilon=0,*parameters):  
    # Extract the wave speed c and the taxis parameter k1
    c = parameters[0]
    k1 = parameters[1]    
    # Return the dynamics of the logistic symmetry
    return array([ -exp(((c-k1)*X[0])/(5))*((5)/(2*(c-k1))),
                  exp(((c-k1)*X[0])/(5))*X[1],
                exp(((c-k1)*X[0])/(5))*(((c-k1)/(5))*X[1]+(3/2)*X[2])])
# Function 5: dX_dt_generalised
# The second order generalised logistic model as a system of two ODEs in terms
# of the pseudo state v=u'. 
def dX_dt_generalised(X, t=0,*parameters):
    # Extract the wave speed c, the first taxis parameter k1
    # and the second taxis parameter k2
    c = parameters[0]
    k1 = parameters[1]
    k2 = parameters[2]    
    # Return the dynamics of the logistic model
    return array([ X[1],
                   -(c-(k1+k2*X[0]))*X[1]-X[0]*((1-X[0])**2)])
# Function 6: Gamma_generalised
# This function plots the symmetry of the logistic model. The first
# Component of this ODE system is xi, the second is eta and the
# third one is eta^{(1)}.
def Gamma_generalised(X, epsilon=0,*parameters):  
    # Extract the wave speed c and the taxis parameter k1
    c = parameters[0]
    k1 = parameters[1]    
    # Return the dynamics of the logistic symmetry
    return array([ -exp(((c-k1)*X[0])/(3))*((3)/(c-k1)),
                  exp(((c-k1)*X[0])/(3))*X[1],
                exp(((c-k1)*X[0])/(3))*(((c-k1)/(3))*X[1]+2*X[2])])
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
# Overall parameters
#=================================================================================
#=================================================================================
# Define the constants k1 and k2 for the taxis
k1 = 1
k2 = 2*sqrt(2)
print("k1=%0.3f"%(k1))
print("k2=%0.3f"%(k2))
# Define the travelling wave speed in terms of these parameters
c_exp = k1 + 5/sqrt(6)
c_log = k1 + 5/sqrt(6)
c_gen = k1 + 3/sqrt(2)
print("c_exp=%0.3f"%(c_exp))
print("c_log=%0.3f"%(c_log))
print("c_gen=%0.3f"%(c_gen))

# The transformation parameter with which the symmetries transform the original
# solution once
epsilon = 0.5
epsilon_vec = linspace(0,epsilon,50)
# Define the initial conditions of the original solution
X_0 = array([0.5, 0])
# Define the travelling wave vector
z_vec = linspace(-1,10,300)
# Define the magical indices which we transform
#magical_indices = arange(0,len(z_vec),round(len(z_vec)/6))
#magical_indices = array([0,5,10,15])
magical_indices = arange(0,140,15)
#=================================================================================
#=================================================================================
# Plotting the symmetries of the exponential model
#=================================================================================
#=================================================================================
# Solve the ODE for the travelling wave equation
X_exp_1, infodict = integrate.odeint(dX_dt_exponential, X_0, z_vec, args=((c_exp,)),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the solution to our lovely system
u_exp, v_exp = X_exp_1.T
# Allocate memory for the transformed states
z_trans_exp = []
u_trans_exp = []
v_trans_exp = []
# Loop over the magical indices
for magical_index in magical_indices:
    # Solve the ODE for the symmetry
    Gamma_exp, infodict = integrate.odeint(Gamma_exponential, array([z_vec[magical_index], u_exp[magical_index], v_exp[magical_index]]), epsilon_vec,full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract the solution to our lovely system
    z_temp, u_temp, v_temp = Gamma_exp.T
    # Save our solutions here
    z_trans_exp.append(z_temp)
    u_trans_exp.append(u_temp)
    v_trans_exp.append(v_temp)    
# Define a transformed z vector
z_hat = linspace(z_trans_exp[0][-1],z_vec[-1],300)
# Plot the transformed solution
X_exp_2, infodict = integrate.odeint(dX_dt_exponential, array([u_trans_exp[0][-1], v_trans_exp[0][-1]]), z_hat, args=((c_exp,)),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the solution to our lovely system
u_hat_exp, v_hat_exp = X_exp_2.T
#=================================================================================
#=================================================================================
# Plotting the symmetries of the logistic model
#=================================================================================
#=================================================================================
# Solve the ODE for the travelling wave equation
X_log_1, infodict = integrate.odeint(dX_dt_logistic, X_0, z_vec, args=((c_log,k1)),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the solution to our lovely system
u_log, v_log = X_log_1.T
# Allocate memory for the transformed states
z_trans_log = []
u_trans_log = []
v_trans_log = []
# Loop over the magical indices
for magical_index in magical_indices:
    # Solve the ODE for the symmetry
    Gamma_log, infodict = integrate.odeint(Gamma_logistic, array([z_vec[magical_index], u_log[magical_index], v_log[magical_index]]), epsilon_vec,args=((c_log,k1)),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract the solution to our lovely system
    z_temp, u_temp, v_temp = Gamma_log.T
    # Save our solutions here
    z_trans_log.append(z_temp)
    u_trans_log.append(u_temp)
    v_trans_log.append(v_temp)
# Define a transformed z vector
z_hat_log = linspace(z_trans_log[0][-1],z_vec[-1],300)
# Plot the transformed solution
X_log_2, infodict = integrate.odeint(dX_dt_logistic, array([u_trans_log[0][-1], v_trans_log[0][-1]]), z_hat_log, args=((c_log,k1)),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the solution to our lovely system
u_hat_log, v_hat_log = X_log_2.T
#=================================================================================
#=================================================================================
# Plotting the symmetries of the generalised logistic model
#=================================================================================
#=================================================================================
# Solve the ODE for the travelling wave equation
X_gen_1, infodict = integrate.odeint(dX_dt_generalised, X_0, z_vec, args=((c_gen,k1,k2)),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the solution to our lovely system
u_gen, v_gen = X_gen_1.T
# Allocate memory for the transformed states
z_trans_gen = []
u_trans_gen = []
v_trans_gen = []
# Loop over the magical indices
for magical_index in magical_indices:
    # Solve the ODE for the symmetry
    Gamma_gen, infodict = integrate.odeint(Gamma_generalised, array([z_vec[magical_index], u_gen[magical_index], v_gen[magical_index]]), epsilon_vec,args=((c_gen,k1)),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract the solution to our lovely system
    z_temp, u_temp, v_temp = Gamma_gen.T
    # Save our solutions here
    z_trans_gen.append(z_temp)
    u_trans_gen.append(u_temp)
    v_trans_gen.append(v_temp)
# Define a transformed z vector
z_hat_gen = linspace(z_trans_gen[0][-1],z_vec[-1],300)
# Plot the transformed solution
X_gen_2, infodict = integrate.odeint(dX_dt_generalised, array([u_trans_gen[0][-1], v_trans_gen[0][-1]]), z_hat_gen, args=((c_gen,k1,k2)),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the solution to our lovely system
u_hat_gen, v_hat_gen = X_gen_2.T    
#=================================================================================
#=================================================================================
# Plot all of our symmetries
#=================================================================================
#=================================================================================
#Define the first figure
f_sym_tv, ax_sym_tv = plt.subplots(1, 3, constrained_layout=True, figsize=(20, 8))
#---------------------------------------------------------------------------------
# Subplot 1 of 3: The exponential model
# The original solution
ax_sym_tv[0].plot(z_vec,u_exp,color=(0/256,68/256,27/256),label="$u(z)$",linewidth=1.0)
# The symmetry
for index in range(len(magical_indices)):
    if index == 0:
        ax_sym_tv[0].plot(z_trans_exp[index], u_trans_exp[index], '--', label="$\\Gamma_{\\epsilon}^{\mathrm{Exponential}}$",color=(0/256,0/256,0/256),linewidth=1.0)
    else:
        ax_sym_tv[0].plot(z_trans_exp[index], u_trans_exp[index], '--',color=(0/256,0/256,0/256),linewidth=1.0)
# The transformed solution        
ax_sym_tv[0].plot(z_hat,u_hat_exp,color=(153/256,216/256,201/256),label="$\\hat{u}(z)$",linewidth=1.0)
# Set a grid and define a legend
ax_sym_tv[0].grid()
ax_sym_tv[0].legend(loc='best',prop={"size":15})
# Set the x-labels and y-labels
ax_sym_tv[0].set_xlabel(xlabel="Travelling wave variable, $z=x-ct$",fontsize=15)
ax_sym_tv[0].set_ylabel(ylabel="Population density, $u(z)$",fontsize=15)
ax_sym_tv[0].set_title(label="Exponential model",fontsize=20)
#---------------------------------------------------------------------------------
# Subplot 2 of 3: The logistic model
# The original solution
ax_sym_tv[1].plot(z_vec,u_log,color=(0/256,68/256,27/256),label="$u(z)$",linewidth=1.0)
# The symmetry
for index in range(len(magical_indices)):
    if index == 0:
        ax_sym_tv[1].plot(z_trans_log[index], u_trans_log[index], '--', label="$\\Gamma_{\\epsilon}^{\mathrm{Logistic}}$",color=(0/256,0/256,0/256),linewidth=1.0)
    else:
        ax_sym_tv[1].plot(z_trans_log[index], u_trans_log[index], '--',color=(0/256,0/256,0/256),linewidth=1.0)
# The transformed solution        
ax_sym_tv[1].plot(z_hat_log,u_hat_log,color=(153/256,216/256,201/256),label="$\\hat{u}(z)$",linewidth=1.0)
# Set a grid and define a legend
ax_sym_tv[1].grid()
ax_sym_tv[1].legend(loc='best',prop={"size":15})
# Set the x-labels and y-labels
ax_sym_tv[1].set_xlabel(xlabel="Travelling wave variable, $z=x-ct$",fontsize=15)
ax_sym_tv[1].set_ylabel(ylabel="Population density, $u(z)$",fontsize=15)
ax_sym_tv[1].set_title(label="Logistic model",fontsize=20)
#---------------------------------------------------------------------------------
# Subplot 3 of 3: The generalised logistic model
# The original solution
ax_sym_tv[2].plot(z_vec,u_gen,color=(0/256,68/256,27/256),label="$u(z)$",linewidth=1.0)
# The symmetry
for index in range(len(magical_indices)):
    if index == 0:
        ax_sym_tv[2].plot(z_trans_gen[index], u_trans_gen[index], '--', label="$\\Gamma_{\\epsilon}^{\mathrm{Generalised}}$",color=(0/256,0/256,0/256),linewidth=1.0)
    else:
        ax_sym_tv[2].plot(z_trans_gen[index], u_trans_gen[index], '--',color=(0/256,0/256,0/256),linewidth=1.0)
# The transformed solution        
ax_sym_tv[2].plot(z_hat_gen,u_hat_gen,color=(153/256,216/256,201/256),label="$\\hat{u}(z)$",linewidth=1.0)
# Set a grid and define a legend
ax_sym_tv[2].grid()
ax_sym_tv[2].legend(loc='best',prop={"size":15})
# Set the x-labels and y-labels
ax_sym_tv[2].set_xlabel(xlabel="Travelling wave variable, $z=x-ct$",fontsize=15)
ax_sym_tv[2].set_ylabel(ylabel="Population density, $u(z)$",fontsize=15)
ax_sym_tv[2].set_title(label="Generalised logistic model",fontsize=20)
#---------------------------------------------------------------------------------
# Show the figure
plt.show()
# Title and saving the figure
f_sym_tv.suptitle('Travelling wave symmetries of models of collective cell migration',fontsize=30,weight='bold');
f_sym_tv.savefig('../Figures/travelling_wave_symmetries.png')

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# PLOTTING THESE TRAVELLING WAVE MODELS IN LATEX
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# THE EXPONENTIAL MODEL
plot_LaTeX_2D(z_vec,u_exp,"../Figures/travelling_waves/Input/exponential.tex","color=exp_1,line width=1pt,","$u(z)$")
plot_LaTeX_2D(z_hat,u_hat_exp,"../Figures/travelling_waves/Input/exponential.tex","color=exp_2,line width=1pt,","$\\hat{u}(z)$")
# The symmetry
for index in range(len(magical_indices)):
    if index == 0:
        plot_LaTeX_2D(z_trans_exp[index], u_trans_exp[index],"../Figures/travelling_waves/Input/exponential.tex","color=black,->,>=latex,densely dashed","$\\Gamma_{\\epsilon}^{\mathrm{Exponential}}$")
    else:
        plot_LaTeX_2D(z_trans_exp[index], u_trans_exp[index],"../Figures/travelling_waves/Input/exponential.tex","color=black,->,>=latex,densely dashed",[])
# THE LOGISTIC MODEL
plot_LaTeX_2D(z_vec,u_log,"../Figures/travelling_waves/Input/logistic.tex","color=log_1,line width=1pt,","$u(z)$")
plot_LaTeX_2D(z_hat_log,u_hat_log,"../Figures/travelling_waves/Input/logistic.tex","color=log_2,line width=1pt,","$\\hat{u}(z)$")
# The symmetry
for index in range(len(magical_indices)):
    if index == 0:
        plot_LaTeX_2D(z_trans_log[index], u_trans_log[index],"../Figures/travelling_waves/Input/logistic.tex","color=black,->,>=latex,densely dashed","$\\Gamma_{\\epsilon}^{\mathrm{Logistic}}$")
    else:
        plot_LaTeX_2D(z_trans_log[index], u_trans_log[index],"../Figures/travelling_waves/Input/logistic.tex","color=black,->,>=latex,densely dashed",[])
# THE GENERALISED LOGISTIC MODEL
plot_LaTeX_2D(z_vec,u_gen,"../Figures/travelling_waves/Input/generalised.tex","color=gen_1,line width=1pt,","$u(z)$")
plot_LaTeX_2D(z_hat_gen,u_hat_gen,"../Figures/travelling_waves/Input/generalised.tex","color=gen_2,line width=1pt,","$\\hat{u}(z)$")
# The symmetry
for index in range(len(magical_indices)):
    if index == 0:
        plot_LaTeX_2D(z_trans_gen[index], u_trans_gen[index],"../Figures/travelling_waves/Input/generalised.tex","color=black,->,>=latex,densely dashed","$\\Gamma_{\\epsilon}^{\mathrm{Generalised}}$")
    else:
        plot_LaTeX_2D(z_trans_gen[index], u_trans_gen[index],"../Figures/travelling_waves/Input/generalised.tex","color=black,->,>=latex,densely dashed",[])

        
# plot_LaTeX_2D(t_lin,u_lin,"../Figures/LaTeX_figures/lift_linear_model/Input/scaling_" + str(case_number) + ".tex","color=u,line width=2pt,","$u(t)$")
# plot_LaTeX_2D(t_lin_scale,u_lin_scale,"../Figures/LaTeX_figures/lift_linear_model/Input/scaling_" + str(case_number) + ".tex","color=u_hat,line width=2pt,","$\\hat{u}(\\hat{t})$")
# plot_LaTeX_2D(t_lin,v_lin,"../Figures/LaTeX_figures/lift_linear_model/Input/scaling_" + str(case_number) + ".tex","color=v,line width=2pt,","$v(t)$")
# plot_LaTeX_2D(t_lin_scale,v_lin_scale,"../Figures/LaTeX_figures/lift_linear_model/Input/scaling_" + str(case_number) + ".tex","color=v_hat,line width=2pt,","$\\hat{v}(\\hat{t})$")
# # Plot the action of the symmetry
# for index in range(0,len(Gamma_S_t)-1):
#     if index == 0:
#         plot_LaTeX_2D(Gamma_S_t[index], Gamma_S_u[index],"../Figures/LaTeX_figures/lift_linear_model/Input/scaling_" + str(case_number) + ".tex","color=black,->,>=latex,densely dashed","$\\Gamma_{2}^{S}$")
#         plot_LaTeX_2D(Gamma_S_t[index], Gamma_S_v[index],"../Figures/LaTeX_figures/lift_linear_model/Input/scaling_" + str(case_number) + ".tex","color=black,->,>=latex,densely dashed",[])        
#     else:
#         plot_LaTeX_2D(Gamma_S_t[index], Gamma_S_u[index],"../Figures/LaTeX_figures/lift_linear_model/Input/scaling_" + str(case_number) + ".tex","color=black,->,>=latex,densely dashed",[])
#         plot_LaTeX_2D(Gamma_S_t[index], Gamma_S_v[index],"../Figures/LaTeX_figures/lift_linear_model/Input/scaling_" + str(case_number) + ".tex","color=black,->,>=latex,densely dashed",[])


