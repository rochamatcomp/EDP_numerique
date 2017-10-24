// Utils to numerical methods for unsteady partial differential equations (PDE)
// Author: Andre ALMEIDA ROCHA
// Professor: M. Bruno Despres
// Instituition: Universite Pierre et Marie Curie

// Update the solution for the periodicity
function [solution] = periodicity(solution, lastSpatial, time)
    solution(lastSpatial + 2, time) = solution(2, time)
    solution(1, time) = solution(lastSpatial + 1, time)
endfunction

// Analytical solver
function [solution] = analyticalSolver(spatials, times, spatialVector, timeVector, speed)
    solution = zeros(spatials, times)
    
    for time = 1:times
        for spatial = 1:spatials
            solution(spatial, time) = initialCondition(spatialVector(spatial) - speed * timeVector(time))
        end
    end
endfunction

// Initialisation of the discrete solution
// (main unknowns from 2 to J+1 + periodicity (2 values: 1 and J+2)
function [solution] = discreteInitialisation(spatials, times)
    solution = zeros(spatials + 2, times)

    // Initialisation to initial time = 1
    time = 1
    for spatial = 1:spatials
        solution(spatial + 1, time) = initialCondition(spatialVector(spatial));
    end
endfunction

// Discrete solver
function [solution]= discreteSolver(scheme, solution, spatials, times, CFL)
    for time = 1:times-1
        solution = periodicity(solution, spatials, time)
    
        for spatial = 2:spatials + 1                
            prev = solution(spatial - 1, time)
            current = solution(spatial, time)
            next = solution(spatial + 1, time)
            
            solution(spatial , time + 1) = scheme(prev, current, next, CFL)
        end
    end
    
    // End periodicity
    discreteSolution = periodicity(solution, spatials, times)
endfunction

// Lax-Wendroff scheme
function [element] = laxWendroff(prev, current, next, CFL)
    element = current - CFL * (next - prev)/2 - CFL**2 * (2*current - prev - next)/2
endfunction

// Centered explicit scheme
function [element] = centeredExplicit(prev, current, next, CFL)
    element = current  - CFL * (next - prev)/2
endfunction

// Implicit scheme
function [element] = implicitScheme(prev, current, next, CFL)
    dg = ones(1:J);
    dgs = ones(1:J-1)*CFL/2;
    Mat=diag(dg,0)-diag(dgs,-1)+diag(dgs,+1);
    Mat(1,J) = -CFL/2 ; Mat(J,1) = CFL/2 
    Mat = sparse(Mat)
    U(2:J+1,n+1)= linsolve(Mat,-U(2:J+1,n))
endfunction
