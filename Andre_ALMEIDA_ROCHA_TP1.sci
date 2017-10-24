//cindy.guichard@upmc.fr

// Numerical methods for unsteady partial differential equations (PDE)
// Author: Andre ALMEIDA ROCHA <andre.almeida_rocha@etu.upmc.fr>
// Professor: M. Bruno Despres
// Instituition: Universite Pierre et Marie Curie

// TP1: Transport equation with constant coefficients
// Space domain: (0, 1). Periodic boundary conditions: u(x+1,t) = u(x,t), period 1

// Includes the definition of utils functions
exec("utils.sci")

// Defines the initial condition (time = 0)
deff('y = initialCondition(x)','y = sin(2*%pi*x)')

// Space mesh size: J in N*
//mesh = 100
mesh = 2

// Speed
speed = 1

// Courant–Friedrichs–Lewy (CFL) condition
cfl = 0.1

// Final time
//finalTime = 0.25
finalTime = 0.05

// Total of columns (spatials)
spatials = mesh

// Spatial step: delta x = 1/J
spatialStep = 1/mesh

// xj = j * delta x (spatial elements)
// Vector of spatial elements
spatialVector = [0:mesh] * spatialStep

// Time step: delta t = cfl * delta x / speed
timeStep = cfl * spatialStep / speed

// tn = n * delta t (time elements)
// Vector of time elements
timeVector = 0:timeStep:finalTime

// Total of lines (times)
times = length(timeVector)

// Periodic solution: u(j,n) = u(j+J,n), j in [0,J], n in N
discreteSolution = discreteInitialisation(spatials, times)
discreteSolution = discreteSolver(laxWendroff, discreteSolution, spatials, times, cfl)

analyticalSolution = analyticalSolver(spatials, times, spatialVector, timeVector, speed)

// The LaTeX representation of the matrix
//str = prettyprint(discreteSolution)
//xstring(0.5, 0.5, str) // Show the representation in a graphic Windows

print(%io(2), discreteSolution)
print(%io(2), analyticalSolution)

exit()
