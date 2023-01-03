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
import matplotlib.pyplot as plt
import numpy as np
# Make the plotting fancy with LaTeX
plt.rcParams['text.usetex'] = True

#=================================================================================
#=================================================================================
# Create a mesh of the unit square
#=================================================================================
#=================================================================================
num_nodes = 1000 # Number of nodes
mesh = UnitSquareMesh(num_nodes, num_nodes, "left") # Generate the mesh
gdim = mesh.geometry().dim() # Get the dimension of the mesh
#=================================================================================
#=================================================================================
# Plot the FEM mesh
#=================================================================================
#=================================================================================
#Define the first figure
fig_FEM_mesh = plt.gcf() # get current figure
fig_FEM_mesh.set_size_inches(20, 8)
#---------------------------------------------------------------------------------
# Subplot 1 of 2: The FEM mesh
# The original solution
plot(mesh)
# Set a grid and define a legend
plt.grid()
# changing the fontsize of yticks
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
# Set the x-labels and y-labels
plt.xlabel("$x$",fontsize=30)
plt.ylabel("$y$",fontsize=30)
plt.title("Mesh of the unit square",fontsize=40,weight='bold')
#---------------------------------------------------------------------------------
# Save the figure
#plt.savefig('../Figures/mesh_unit_square.png',dpi = 100)
# Show the figure
plt.show()
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
