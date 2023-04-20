#=================================================================================
#=================================================================================
# Script:"wound_healing_FEM_simulation"
# Date: 2023-01-03
# Implemented by: Johannes Borgqvist
# Description:
# We try to implement a FEM simulation of wound healing using the model we
# have constructed from a symmetry.
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
# Most good things are contained in FEniCS
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate # For solving ODEs.
from mpl_toolkits.axes_grid1 import make_axes_locatable
# Make the plotting fancy with LaTeX
plt.rcParams['text.usetex'] = True
#=================================================================================
#=================================================================================
# Functions
#=================================================================================
#=================================================================================
# Function 1: dX_dt_ICs
# This function defines the ODE system on which we base the initial conditions
# of the pseudo state v=u'. 
def dX_dt_ICs(X, t=0,*parameters):   
    # Define our three terms
    diffusion = X[1]**2/X[0]
    taxis = (-2*X[0]+5)*X[1]
    cell_growth = X[0]*(1-X[0])*(X[0]+3)
    # Return the dynamics of the exponential model
    return np.array([ X[1],
                   -(diffusion+taxis+cell_growth)])
# Function 2: cell_migration_IC
# This function maps a travelling wave solution onto the unit square in
# order to create an initial condition (IC) for the continumm PDE
def cell_migration_IC(c_x,c_y,z_vec,u_z,H,u_prev):
    # Create a dofmap
    dofmap = H.dofmap()
    # We want to save all of the dofs as well
    dofs = []
    # Loop over the cells in the mesh and add the dofs in the respective list
    for cell in cells(mesh):
        # Add all the dofs in the sphere
        dofs.extend(dofmap.cell_dofs(cell.index()))            
    # Find the unique dofs
    dofs = list(set(dofs))
    # Get the actual spatial coordinates as well
    coordinates = H.tabulate_dof_coordinates()
    # Loop through the dofs and the coordinates and save the initial conditions
    for dof, coordinate in zip(dofs, coordinates):
        # Calculate the z-value
        z_val = c_x*coordinate[0]+c_y*coordinate[1]
        # Calculate the corresponding u-value which we save at the correct
        # position in the mesh
        u_prev.vector()[dof] = np.interp(z_val, z_vec, u_z)



#=================================================================================
#=================================================================================
# Create a mesh of the unit square
#=================================================================================
#=================================================================================
# Define the initial conditions of the original solution
X_0 = np.array([0.80, 0])
# Define the travelling wave vector
z_vec = np.linspace(-1,5,300)
# Solve the ODE for the travelling wave equation
X_IC, infodict = integrate.odeint(dX_dt_ICs, X_0, z_vec,full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the solution to our lovely system
u_z, v_z = X_IC.T
#=================================================================================
#=================================================================================
# Plot the travelling wave solution
#=================================================================================
#=================================================================================
#Define the first figure
fig_TW = plt.figure(1) # get current figure
fig_TW.set_size_inches(10, 10)
#---------------------------------------------------------------------------------
# The travelling wave solution
plt.plot(z_vec,u_z)
# Set a grid and define a legend
plt.grid()
# changing the fontsize of yticks
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
# Set the x-labels and y-labels
plt.xlabel("Travelling wave variable, $z$",fontsize=30)
plt.ylabel("Population density, $u(z)$",fontsize=30)
plt.title("Travelling wave solution with density dependent diffusion",fontsize=40,weight='bold')
#plt.show()


#=================================================================================
#=================================================================================
# Create a mesh of the unit square
#=================================================================================
#=================================================================================
num_nodes = 25 # Number of nodes
mesh = UnitSquareMesh(num_nodes, num_nodes, "left") # Generate the mesh
#=================================================================================
#=================================================================================
# Illustrate the initial conditions for a certain choice of wave speeds 
#=================================================================================
#=================================================================================
# Collect everything in the function space (name it H for a Hilbert space)
H = FunctionSpace(mesh, 'P', 2)
# Define our wave speed variables
c_x = 1/np.sqrt(2)
c_y = 1/np.sqrt(2)
c_t = 1
# Define a function in which we will save the initial condition
u_prev = Function(H) # Previous time step in the FD time stepping scheme
# Calculate the initial conditions
cell_migration_IC(c_x,c_y,z_vec,u_z,H,u_prev)
# Set the time to 0
t = 0
#=================================================================================
#=================================================================================
# Find the wave coordinates
#=================================================================================
#=================================================================================
# Find the medium value of the wave
wave_front = 0.65*max(u_z)
# Find the z-coordinate of the wave front
z_front = np.interp(wave_front, np.flip(u_z), np.flip(z_vec))
# Now, we make a scatter vector where we plot x, vs y of the front
x_front = np.linspace(0,1)
y_front = np.array([((z_front-c_x*x+c_t*t)/(c_y)) for x in x_front])
#=================================================================================
#=================================================================================
# Plot the initial conditions
#=================================================================================
#=================================================================================
#Define the second figure
fig_IC = plt.figure(3) # get current figure
fig_IC.set_size_inches(10, 10)
ax = plt.gca()
#---------------------------------------------------------------------------------
# The FEM mesh of the unit square
c=plot(u_prev,mode='color',vmin=0.32,vmax=0.65)
# Plot the wave front as well
plt.plot(x_front,y_front,color=(1,1,1),linewidth=4)
# Set the limits to the unit square
plt.xlim([0,1])
plt.ylim([0,1])
# Set a grid and define a legend
plt.grid()
# changing the fontsize of yticks
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
# Set the x-labels and y-labels
plt.xlabel("$x$",fontsize=40)
plt.ylabel("$y$",fontsize=40)
plt.title("Initial conditions",fontsize=40,weight='bold')
plt.savefig('../Figures/initial_condition.png',dpi = 100)
plt.savefig('../Figures/initial_conditions/initial_condition.png',dpi = 500)
plt.show()







