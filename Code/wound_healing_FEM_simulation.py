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
from dolfin import *
import matplotlib.pyplot as plt
# Make the plotting fancy with LaTeX
plt.rcParams['text.usetex'] = True





mesh = UnitSquareMesh(10, 10, "left")
print("Plotting a UnitIntervalMesh")



#=================================================================================
#=================================================================================
# Plot all of our symmetries
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
plt.savefig('../Figures/mesh_unit_square.png',dpi = 100)
# Show the figure
plt.show()

