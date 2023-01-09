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
# Function 3: solve_PDE_cell_migration
# This function solves the continuum PDE modelling cell migration
def solve_PDE_cell_migration(num_nodes,mesh,c_x,c_y,z_vec,u_z,parameters):
    #---------------------------------------------------------------------------------
    # Step 0 out of 4: Extract parameters
    #---------------------------------------------------------------------------------
    T = parameters[0] # End time for simulation
    v_min = parameters[1] # Lower limit for colour bar
    v_max = parameters[2] # Upper limit for colour bar
    T1 = parameters[3] # First time point we save the figure of
    T2 = parameters[4] # Second time point we save the figure of
    T3 = parameters[5] # Third time point we save the figure of
    #---------------------------------------------------------------------------------
    # Step 1 out of 4: Setup function spaces
    #---------------------------------------------------------------------------------
    # Collect everything in the function space (name it H for a Hilbert space)
    H = FunctionSpace(mesh, 'P', 2)
    # Define our functions needed and testfunctions in the FD-FEM time stepping scheme    
    u_prev = Function(H) # Previous time step in the FD time stepping scheme
    # Define our analytical solution
    u = Function(H)
    # Define our lovely test function
    v = TestFunction(H)    
    #---------------------------------------------------------------------------------
    # Step 2 out of 4: Setup initial conditions
    #---------------------------------------------------------------------------------
    # Calculate the initial conditions
    cell_migration_IC(c_x,c_y,z_vec,u_z,H,u_prev)
    # Save the initial condition and look at it in ParaView
    vtkfile_u = File("../Output/cell_migration_" + "c_x_" + str(round(c_x,3)).replace(".","p") + "_c_y_" + str(round(c_y,3)).replace(".","p") + "/u.pvd")
    t = 0.0 # The time is zero to begin with
    # WE ALSO SAVE THE VERY LAST ITERATION WHEN ALL THE TIME STEPPING IS DONE.
    u_prev.rename("Population density, $u(\mathbf{x},t)$","u")
    vtkfile_u << (u_prev, t)
    # We want to plot the initial conditions as well
    #---------------------------------------------------------------------------------
    #Define the second figure
    fig_temp = plt.figure() # get current figure
    fig_temp.set_size_inches(10, 10)
    ax = plt.gca()
    # The FEM mesh of the unit square
    c=plot(u_prev,mode='color',vmin=v_min,vmax=v_max)
    # Set a grid and define a legend
    plt.grid()
    # changing the fontsize of yticks
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    # Set the x-labels and y-labels
    plt.xlabel("$x$",fontsize=40)
    plt.ylabel("$y$",fontsize=40)
    plt.title("$t=" + str(round(t,2)) + "$",fontsize=40,weight='bold')
    plt.savefig("../Figures/PDE_simulations/Input/cx_" + str(round(c_x,3)).replace(".","p") + "_" + str(round(c_y,3)).replace(".","p") + "_" + str(round(t,2)).replace(".","p") + ".png",dpi = 500)
    #---------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------
    # Step 3 out of 4: Setup numerical parameters
    #---------------------------------------------------------------------------------
    # Define numerical parameters for the FD-scheme in time
    delta_x = 1/num_nodes # The step length in the spatial dimension
    dt = ((delta_x**2)/(8)) # Time step in finite difference
    # Re-define this as a constant for the FEM solver    
    k = Constant(dt) # For the fem solver as well
    # Define an iterator for the time stepping keeping track of
    # how many iterations that has passed
    t_it = 0
    # Also, define the current time outside the loop, so that we can save
    # the final concentration profile after the looping is done.
    t = 0
    # Previous time step
    t_prev = 0
    # Save the concentration profile when the time has increased with a value 0.5
    save_iteration_vtk_original = T/60
    save_iteration_vtk = save_iteration_vtk_original
    #--------------------------------------------------------------------------------
    # Step 4 out of 4: FD time stepping, a simple Euler forward, in time FEM in space
    #--------------------------------------------------------------------------------
    # Surpress annoying printing from the FEM solver
    set_log_active(False)
    # Start the time stepping
    while t < T:
        #--------------------------------------------------------------------------------
        # Update the time stepping details
        #--------------------------------------------------------------------------------    
        # Save the previous time step
        t_prev = t
        # Update current time and iteration number 
        t_it += 1
        t += dt
        #--------------------------------------------------------------------------------
        # FEM solution for current time step
        #--------------------------------------------------------------------------------
        # Define the individual terms in the VF
        Accumulation = (u-u_prev)*v*dx
        Diffusion = u*dot(grad(u), grad(v))*dx
        Taxis = ((-(2*u**2)+(5*u)-1)/(c_x+c_y))*(grad(u)[0]+grad(u)[1])*v*dx
        Cell_growth = (u**2)*(1-u)*(u+3)*v*dx
        # Define our function that we want to solve for
        F = Accumulation + k*(Diffusion-Taxis-Cell_growth)
        # Solve for u
        solve(F==0,u)
        #--------------------------------------------------------------------------------
        # Save and update solution
        #--------------------------------------------------------------------------------    
        # Save and check the solution (every whatever iteration)
        if  (t_prev < save_iteration_vtk) and (t > save_iteration_vtk):
            # Prompt to the user
            print("\t\tIteration %d, t\t=\t%0.3f out of %0.3f"%(t_it,t,T))    
            # Increase the iterations
            save_iteration_vtk += save_iteration_vtk_original
            # Re--name the title on the colour bar
            u.rename("Population density, $u(\mathbf{x},t)$","u")
            # Save the components in the data files
            vtkfile_u << (u, t)
        # Update old solution
        u_prev.assign(u)
        # We also want to save some nice figures here
        if ((t_prev < T1) and (t > T1)) or ((t_prev < T2) and (t > T2)) or ((t_prev < T3) and (t > T3)):
            #Define the second figure
            fig_temp = plt.figure() # get current figure
            fig_temp.set_size_inches(10, 10)
            ax = plt.gca()
            #---------------------------------------------------------------------------------
            # The FEM mesh of the unit square
            c=plot(u,mode='color',vmin=v_min,vmax=v_max)
            # Set a grid and define a legend
            plt.grid()
            # changing the fontsize of yticks
            plt.xticks(fontsize=30)
            plt.yticks(fontsize=30)
            # Set the x-labels and y-labels
            plt.xlabel("$x$",fontsize=40)
            plt.ylabel("$y$",fontsize=40)
            plt.title("$t=" + str(round(t,2)) + "$",fontsize=40,weight='bold')
            plt.savefig("../Figures/PDE_simulations/Input/cx_" + str(round(c_x,3)).replace(".","p") + "_" + str(round(c_y,3)).replace(".","p") + "_" + str(round(t,2)).replace(".","p") + ".png",dpi = 500)          
    #--------------------------------------------------------------------------------
    # Save solution when t=T
    #--------------------------------------------------------------------------------
    # Re--name the title on the colour bar
    u.rename("Population density, $u(\mathbf{x},t)$","u")
    # Save the components in the data files
    vtkfile_u << (u, t)
    #---------------------------------------------------------------------------------
    #Define the second figure
    fig_temp = plt.figure() # get current figure
    fig_temp.set_size_inches(10, 10)
    ax = plt.gca()
    # The FEM mesh of the unit square
    c=plot(u,mode='color',vmin=v_min,vmax=v_max)
    # Set a grid and define a legend
    plt.grid()
    # changing the fontsize of yticks
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    # Set the x-labels and y-labels
    plt.xlabel("$x$",fontsize=40)
    plt.ylabel("$y$",fontsize=40)
    plt.title("$t=" + str(round(t,2)) + "$",fontsize=40,weight='bold')
    plt.savefig("../Figures/PDE_simulations/Input/cx_" + str(round(c_x,3)).replace(".","p") + "_" + str(round(c_y,3)).replace(".","p") + "_" + str(round(t,2)).replace(".","p") + ".png",dpi = 500)
    #---------------------------------------------------------------------------------    
    # Prompt to the user
    print("\n\tSimulations are finished!\n\n")        
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
plt.show()
plot_LaTeX_2D(z_vec,u_z,"../Figures/initial_conditions/Input/u.tex","color=exp_1,line width=2pt,",[])

#=================================================================================
#=================================================================================
# Create a mesh of the unit square
#=================================================================================
#=================================================================================
num_nodes = 25 # Number of nodes
mesh = UnitSquareMesh(num_nodes, num_nodes, "left") # Generate the mesh
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
# Illustrate the initial conditions for a certain choice of wave speeds 
#=================================================================================
#=================================================================================
# Collect everything in the function space (name it H for a Hilbert space)
H = FunctionSpace(mesh, 'P', 2)
# Define our wave speed variables
c_x = 1/np.sqrt(2)
c_y = 1/np.sqrt(2)
# Define a function in which we will save the initial condition
u_prev = Function(H) # Previous time step in the FD time stepping scheme
# Calculate the initial conditions
cell_migration_IC(c_x,c_y,z_vec,u_z,H,u_prev)
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
#=================================================================================
#=================================================================================
# SOLVE THE PDE FOR THE GIVEN PARAMETERS
#=================================================================================
#=================================================================================
# Re-define the mesh with a higher node density
num_nodes = 32 # Number of nodes
mesh = UnitSquareMesh(num_nodes, num_nodes, "left") # Generate the mesh
# Define some hand picked parameters for the visualisation
parameters = [0.20, 0.32, 0.65, 0.05, 0.10, 0.15]
# First solution to PDE problem
solve_PDE_cell_migration(num_nodes,mesh,c_x,c_y,z_vec,u_z,parameters)
# Re-define our wave speeds
c_x = -1/2
c_y = np.sqrt(3)/2
# Define some hand picked parameters for the visualisation
parameters = [0.10, 0.44, 0.75, 0.025, 0.05, 0.075]
# Solve the system of PDEs again with the new wave speed parameters
solve_PDE_cell_migration(num_nodes,mesh,c_x,c_y,z_vec,u_z,parameters)
