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
plot_LaTeX_2D(z_vec,u_z,"../Figures/initial_conditions/Input/u.tex","color=exp_1,line width=2pt,",[])

#=================================================================================
#=================================================================================
# Create a mesh of the unit square
#=================================================================================
#=================================================================================
num_nodes = 25 # Number of nodes
mesh = UnitSquareMesh(num_nodes, num_nodes, "left") # Generate the mesh
gdim = mesh.geometry().dim() # Get the dimension of the mesh
#=================================================================================
#=================================================================================
# Plot the FEM mesh
#=================================================================================
#=================================================================================
#Define the second figure
fig_FEM_mesh = plt.figure(2) # get current figure
fig_FEM_mesh.set_size_inches(10, 10)
#---------------------------------------------------------------------------------
# The FEM mesh of the unit square
plot(mesh)
# Set a grid and define a legend
plt.grid()
# changing the fontsize of yticks
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
# Set the x-labels and y-labels
plt.xlabel("$x$",fontsize=40)
plt.ylabel("$y$",fontsize=40)
plt.title("$" + str(num_nodes) + "\\times" + str(num_nodes) + "$--mesh of the unit square",fontsize=40,weight='bold')
#---------------------------------------------------------------------------------
# Save the figure
plt.savefig('../Figures/mesh_unit_square.png',dpi = 100)
plt.savefig('../Figures/initial_conditions/mesh_unit_square.eps',dpi = 500)
# Show the figure
#plt.show()
#=================================================================================
#=================================================================================
# DEFINE HILBERT SPACE AND SET INITIAL CONDITIONS
#=================================================================================
#=================================================================================
# STEP 1 OUT OF 2: DEFINE THE FUNCTION SPACE
# Define the basis finite elements and their order. We pick linear basis functions. 
#P1 = FiniteElement("P", mesh.ufl_cell(), 1)
# Define a mixed element since we have a coupled system of two PDEs
#element = MixedElement([P1])
# Collect everything in the function space (name it H for a Hilbert space)
H = FunctionSpace(mesh, 'P', 2)
# Define our wave speed variables
c_x = 1/np.sqrt(2)
c_y = 1/np.sqrt(2)
# Define a function in which we will save the initial condition
u_prev = Function(H) # Previous time step in the FD time stepping scheme
# Calculate the initial conditions
cell_migration_IC(c_x,c_y,z_vec,u_z,H,u_prev)
# Save the initial condition and look at it in ParaView
vtkfile_u = File("../Output/wound_healing/u.pvd")
t = 0.0 # The time is zero to begin with
# WE ALSO SAVE THE VERY LAST ITERATION WHEN ALL THE TIME STEPPING IS DONE.
u_prev.rename("Population density, $u(\mathbf{x},t)$","u")
vtkfile_u << (u_prev, t)
#=================================================================================
#=================================================================================
# Setup the variational formulation
#=================================================================================
#=================================================================================
# Define our analytical solution
u = Function(H)
# Define our lovely test function
v = TestFunction(H)
# Define our function that we want to solve for
F = dot(grad(u), grad(v))*dx + (grad(u)[0]+grad(u)[1])*v*dx


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
c=plot(u_prev,mode='color',vmin=0,vmax=1)
# For the colourbar
#cax = fig_IC.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
# The colourbar in the correct size
#cb = plt.colorbar(c,ticks=[0, 0.2, 0.8, 1],cax=cax)
#cb.ax.tick_params(labelsize=30)
#cb.ax.set_ylabel("$u(x,y,t=0)$", rotation=270,loc='center',fontsize=40)
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
