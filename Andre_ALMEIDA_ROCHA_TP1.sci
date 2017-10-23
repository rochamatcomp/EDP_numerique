//cindy.guichard@upmc.fr

// Numerical methods for unsteady partial differential equations (PDE)
// Author: Andre ALMEIDA ROCHA
// Professor: M. Bruno Despres
// Instituition: Universite Pierre et Marie Curie

// TP1: Transport equation with constant coefficients

// Space domain: (0, 1). Periodic boundary conditions: u(x+1,t) = u(x,t), period 1

deff('y = initial_condition(x)','y = sin(2*%pi*x)')

function [x, y]= lax_wendroff(spatials, times, CFL)
    // Lax-Wendroff
    for spatial = 2:spatials + 1
        U(j,n+1)= U(j,n)  -CFL*(U(j+1,n) - U(j-1,n))/2 - CFL**2*(2*U(j,n) -U(j-1,n)-U(j+1,n) )/2
    end
endfunction
    
    

// Space mesh size: J in N*
mesh_size = 100

// Spatial step: delta x = 1/J
spatial_step = 1/mesh_size

// xj = j * delta x (spatial elements)
// Vector of spatial elements
spatial_vector = [0:mesh_size] * spatial_step

// Speed 
speed = 1

// CFL
CFL = 0.1

// Final time
final_time = 0.25

// Time step
time_step = CFL*spatial_step/speed

// tn = n * delta t (time elements)
// Vector of time elements
time_vector = 0:time_step:final_time

// Total of lines (times)
times = length(time_vector)

// Total of columns (spatials)
spatials = mesh_size

// Discrete solution (main unknowns from 2 to J+1 + periodicity (2 values: 1 and J+2)
discrete_solution = zeros(spatials + 2, times)

// Analytical solution
analytical_solution = zeros(spatials, times)

// Initialisation to initial time = 1
time = 1
for spatial = 1:spatials
     discrete_solution(spatial + 1, time) = initial_condition(spatial_vector(spatial));
end

// Periodic solution: u(j,n) = u(j+J,n), j in [0,J], n in N
for time = 1:times-1

    // Periodicity
    //discrete_solution(spatials + 2, time) = discrete_solution(2, time)
    //discrete_solution(1, time) = discrete_solution(spatial + 1, time)
    
    print(%io(2), discrete_solution(spatials + 1, time), discrete_solution(1, time))
end


for time = 1:times
    for spatial = 1:spatials
        analytical_solution(spatial, time) = initial_condition(spatial_vector(spatial) - speed * time_vector(time))        
    end
end

//print(%io(2), solution)
