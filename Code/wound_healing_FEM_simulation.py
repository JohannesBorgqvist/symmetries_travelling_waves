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
#=================================================================================
#=================================================================================
# Create a mesh of the unit square
#=================================================================================
#=================================================================================
# Define the initial conditions of the original solution
X_0 = np.array([0.5, 0])
# Define the travelling wave vector
z_vec = np.linspace(-1,5,300)
# Solve the ODE for the travelling wave equation
X_IC, infodict = integrate.odeint(dX_dt_ICs, X_0, z_vec,full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the solution to our lovely system
u_z, v_z = X_IC.T
#print(u_z)
#=================================================================================
#=================================================================================
# Plot the travelling wave solution
#=================================================================================
#=================================================================================
#Define the first figure
fig_TW = plt.figure(1) # get current figure
fig_TW.set_size_inches(20, 8)
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
gdim = mesh.geometry().dim() # Get the dimension of the mesh
#=================================================================================
#=================================================================================
# Plot the FEM mesh
#=================================================================================
#=================================================================================
#Define the second figure
fig_FEM_mesh = plt.figure(2) # get current figure
fig_FEM_mesh.set_size_inches(20, 8)
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
#plt.savefig('../Figures/mesh_unit_square.png',dpi = 100)
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
# STEP 2 OUT OF 2: SET THE INITIAL CONDITIONS
# Create a dofmap
dofmap = H.dofmap()
# Classify the dofs as either in the hole or on the rest of the sphere
dofs = []
# Loop over the cells in the mesh and add the dofs in the respective list
for cell in cells(mesh):
    # Add all the dofs in the sphere
    dofs.extend(dofmap.cell_dofs(cell.index()))            
# Find the unique dofs
dofs = list(set(dofs))
# Find all the dofs (i.e. the indices of a particular spatial coordinate in the mesh)
#dofs = dofmap.dofs()
#dofs = dofmap.cell_dofs()
# Get the actual coordinates as well
coordinates = H.tabulate_dof_coordinates()
# Define a function which we would like to save the initial condition in
u_prev = Function(H) # Previous time step in the FD time stepping scheme
# Define an elliptic wound
a = 0.4
b = 0.5
# Define a wound threshold
wound_threshold = 0.2
# Loop through the dofs and the coordinates and save the initial conditions
for dof, coordinate in zip(dofs, coordinates):
    # Calculate the ellipse coordinate
    wound_coordinate = ((coordinate[0]-0.5)**2/a**2)+((coordinate[1]-0.5)**2/b**2)
    # If we are outside the ellipse we put it to 1 otherwise 0
    if wound_coordinate>1:
        u_prev.vector()[dof] = 1
    elif wound_coordinate>wound_threshold and wound_coordinate<1.0:
        u_prev.vector()[dof] = exp(1/((1-wound_threshold)**2))*exp(1/((1-wound_coordinate)**2-(1-wound_threshold)**2))
    elif wound_coordinate<wound_threshold:
        u_prev.vector()[dof] = 0
        
# Save the initial condition and look at it in ParaView
vtkfile_u = File("../Output/wound_healing/u.pvd")
t = 0.0 # The time is zero to begin with
# WE ALSO SAVE THE VERY LAST ITERATION WHEN ALL THE TIME STEPPING IS DONE.
u_prev.rename("Population density, $u(\mathbf{x},t)$","u")
vtkfile_u << (u_prev, t)


#=================================================================================
#=================================================================================
# Plot the initial conditions
#=================================================================================
#=================================================================================
#Define the second figure
fig_IC = plt.figure(3) # get current figure
fig_IC.set_size_inches(20, 8)
#fig_IC.yaxis.set_ticks_position('left')
#---------------------------------------------------------------------------------
# The FEM mesh of the unit square
c=plot(u_prev,mode='color')
cb = plt.colorbar(c)
cb.ax.tick_params(labelsize=30)
cb.ax.set_ylabel("$u(x,y,t=0)$", rotation=270,fontsize=40)
# Set a grid and define a legend
plt.grid()
# changing the fontsize of yticks
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
# Set the x-labels and y-labels
plt.xlabel("$x$",fontsize=40)
plt.ylabel("$y$",fontsize=40)
plt.title("Initial conditions",fontsize=40,weight='bold')
plt.show()
