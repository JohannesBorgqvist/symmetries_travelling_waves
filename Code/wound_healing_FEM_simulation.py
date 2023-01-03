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
f_FEM_mesh, ax_FEM_mesh = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
ax_FEM_mesh[0].plot(mesh)
# Set a grid and define a legend
ax_sym_tv[0].grid()
ax_sym_tv[0].legend(loc='best',prop={"size":3})
# Set the x-labels and y-labels
ax_sym_tv[0].set_xlabel(xlabel="$x$",fontsize=5)
ax_sym_tv[0].set_ylabel(ylabel="$y$",fontsize=5)
ax_sym_tv[0].set_title(label="Mesh of the unit square",fontsize=10)
plt.show()

