// Numerical methods for unsteady partial differential equations (PDE)
// Author: Andre ALMEIDA ROCHA <andre.almeida_rocha@etu.upmc.fr>
// Professor: M. Bruno DESPRES
// Instituition: Universite Pierre et Marie Curie - Paris 6

// Homework B004 2017
// 2. Implementation and numerical tests

// a) Scheme

// Gets the flux value.
//
// Inputs 
//  position: 1 if the spatial index of the flux is j + 1/2, -1 if j - 1/2.
//  current: current element of the solution, in current space and fixed time. u(j, n).
//  prev: previous element of the solution, in previous space and fixed time. u(j-1, n).
//  next: next element of the solution, in next space and fixed time. u(j+1, n).
//  speed: advection speed. (a).
//  spatialStep: spatial step of the mesh. (delta x).
//  cflCondition: Courant–Friedrichs–Lewy condition. (a * delta t / delta x).
//  shiftFactor: factor multiplicative of the shift of the derivative. (k).
//
// Outputs
//  flux: flux value.
function [flux] = getFlux(position, current, prev, next, speed, spatialStep, cflCondition, shiftFactor)
    // Shift of the derivatives
    shift = (next - current) / spatialStep
    
    flux = speed * (current + position * 0.5 * (1 - cflCondition) * (current - prev)) - shiftFactor * shift
endfunction

// Calculates the element, in the next time, of the solution by the scheme.
// 
// Inputs
//  current: current element of the solution, in current space and fixed time. u(j, n).
//  prev: previous element of the solution, in previous space and fixed time. u(j-1, n).
//  next: next element of the solution, in next space and fixed time. u(j+1, n).
//  speed: advection speed. (a).
//  spatialStep: spatial step of the mesh. (delta x).
//  timeStep: time step of the mesh. (delta t).
//  shiftFactor: factor multiplicative of the shift of the derivative. (k).
// 
// Outputs
//  element: the element, in the next time, of the solution by the scheme.
function [element] = schemeExercise2a(current, prev, next, speed, spatialStep, timeStep, shiftFactor)
    // Courant–Friedrichs–Lewy condition
    cflCondition = speed * timeStep / spatialStep
    
    fluxNext = getFlux(1, current, prev, next, speed, spatialStep, cflCondition, shiftFactor)
    fluxPrev = getFlux(-1, current, prev, next, speed, spatialStep, cflCondition, shiftFactor)
    element = current + timeStep / spatialStep * (fluxPrev - fluxNext)
endfunction

// b) First test. Consider the following data: a = 2, k = 0, final time = 1.
//    and the delta x, delta t chosen.

// Advection speed. (a).
speed = 2

// Factor multiplicative of the shift of the derivative. (k).
shiftFactor = 0

// Spatial mesh size: J in N*.
spatials = 4

// Spatial step: delta x = 1/J.
spatialStep = 1 / spatials

// Vector of spatial elements (xj = j * delta x).
spatialVector = [0:spatials] * spatialStep

// Final time to the simulation. 
finalTime = 1

// Time step
timeStep = 0.25

//// Courant–Friedrichs–Lewy (CFL) condition
//cflCondition = 0.1
//
//// Time step: delta t = CFL * delta x / speed
//timeStep = cflCondition * spatialStep / speed

// Vector of time elements (tn = n * delta t)
timeVector = 0:timeStep:finalTime

// Time mesh size
times = length(timeVector)

// Redefines the interval as open interval (changes the limits with epsilon)
//epsilon = 0.00001
//spatialVector(1) = spatialVector(1) + epsilon
//spatialVector(spatials+1) = spatialVector(spatials+1) - epsilon

//print(%io(2), spatialVector)
//print(%io(2), timeVector)

solution = zeros(spatials, times)
for time = 1:times
    for space = 1:spatials
        value = spatialVector(space) - speed * timeVector(time)
        solution(space, time) = atan(value) / %pi + 0.5
    end
end

// Update the solution for the periodicity
function [solution] = periodicity(solution, lastSpatial, time)
    solution(lastSpatial + 2, time) = solution(2, time)
    solution(1, time) = solution(lastSpatial + 1, time)
endfunction

//solution = periodicity(solution, spatials, time)


// b1. Implementation of the analytical solution, for first initial condition:
//     u0 = 1, if x dans I = (0.2, 0.3) union (0.6, 0.8)
//          0, if x dans (0, 1) \ I

// Defines the first initial condition (time = 0).
// Inputs
//  spatialElement: value in space position.
// Outputs
//  solution: value of the initial condition by the spatial element.
function [solution] = initialCondition1(spatialElement)
    // Mapping spatial element at the valid domain, interval (0,1).
    // Function $$f: \mathbb{R} \to (0,1), \quad x \mapsto 1/\pi \arctan(x) + 1/2 $$, 
    // as inverse function of the bijection $$ g: (0,1) \to \mathbb{R}, \quad y \mapsto tan(\pi y - \pi/2)$$.
    elementOnDomain = atan(spatialElement) / %pi + 0.5
    
    if (elementOnDomain > 0.2 & elementOnDomain < 0.3) | (elementOnDomain > 0.6 & elementOnDomain < 0.8) then
        solution = 1
    elseif (elementOnDomain > 0 & elementOnDomain <= 0.2) | (elementOnDomain >= 0.3 & elementOnDomain <= 0.6) | (elementOnDomain >= 0.8 & elementOnDomain < 1) then
        solution = 0
    else
        disp(elementOnDomain, "The initial condition is not defined out of the openinterval (0,1).")
        abort
    end
endfunction

// Analytical solver.
//
// Inputs
//  spatials: total of spatials partitions.
//  times: total of times partitions.
//  spatialVector: spatial vector with the spatial elements.
//  timeVector: time vector with the time elements.
//  speed: advection speed. (a).
//
// Outputs
//  solution: analytical solution.
function [solution] = analyticalSolver(spatials, times, spatialVector, timeVector, speed)
    // Initialisation of the solution
    solution = zeros(spatials, times)
    
    for time = 1:times
        for space = 1:spatials
            element = spatialVector(space) - speed * timeVector(time)
            solution(space, time) = initialCondition1(element)
        end
    end
endfunction

analyticalSolution = analyticalSolver(spatials, times, spatialVector, timeVector, speed)
disp(analyticalSolution)

// b2. L-infinite stability of the scheme



// b3. The numerical convergence order of the scheme in L2 norm,
//     for the second initial condition:
//     u0 (x) = alpha sin(2 pi x) + beta cos(2 pi x), alpha, beta dans R.

// Defines the second initial condition (time = 0).
// Inputs
//  spatialElement: value in space position.
// Outputs
//  solution: value of the initial condition by the spatial element. 
function [solution] = initialCondition2(spatialElement, alphaFactor, betaFactor)
    solution = alphaFactor * sin(2 * %pi * spatialElement) + betaFactor * cos(2 * %pi * spatialElement)
endfunction


// c) Second test. Consider following data: a = 0, k = 1, final time = 1,
//    alpha = 1, beta = 0. The analytical solution is:
// u(x, t) = (alpha sin(2 pi (x − at)) + beta cos(2 pi (x − at))) exp(−4 pi^2 t)

// c1. The measurement of the numerical convergence order in L2 norm,
//     for 2 delta t / delta x^2 = 1/2. Comparation with the theory.


// c2. The measurement of the numerical error, for delta t / delta x^2 = 1/6.
//     Numerical evidence that the error is 4th order in space.
//     Explanation of the super-convergence name.
