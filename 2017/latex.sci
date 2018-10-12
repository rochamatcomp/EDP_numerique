txt = gce();
txt.font_size = 5;
txt.font_style = 0;

equation = ["$\overbrace{Scilab}$" "n''est ";"pas" "$\underbrace{Matlab}$"]
xstring(0, 0.5, equation)

spatialElement = 10

disp("Initial condition not definied out the interval (0, 1).", spatialElement)
abort
