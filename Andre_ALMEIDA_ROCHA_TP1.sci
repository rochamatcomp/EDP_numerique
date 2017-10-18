//cindy.guichard@upmc.fr

// Numerical methods for unsteady partial differential equations (PDE)
// Author: Andre ALMEIDA ROCHA
// Professor: M. Bruno Despres
// Instituition: Universite Pierre et Marie Curie

// TP1: Transport equation with constant coefficients

// Space domain: (0, 1). Periodic boundary conditions: u(x+1,t) = u(x,t), period 1

// Space mesh size: J in N*
mesh_size = 100

// Spatial step: delta x = 1/J
spatial_step = 1/mesh_size

// xj = j * delta x (spatial elements)
// Vector of spatial elements
spatials = [0:mesh_size] * spatial_step

// Speed 
speed = 1

// CFL
CFL = 0.1

// time step
time_step = CFL*dx/speed

// tn = n * delta t (time elements)
// Vector of time elements
times = [0:] * step_spatial

print(%io(2), X)


// Periodic solution: u(j,n) = u(j+J,n), j in [0,J], n in N
