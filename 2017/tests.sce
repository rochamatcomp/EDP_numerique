function [solution] = initialCondition2(spatialElement, alpha, beta)
    solution = alpha * sin(2 * %pi * spatialElement) + beta * cos(2 * %pi * spatialElement)
endfunction

solution = initialCondition2(2, 0.5, 0.5)

disp(solution)

exit()
