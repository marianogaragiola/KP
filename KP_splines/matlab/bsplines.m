clear all
close all

pkg load nurbs

num_intervalos   = 5; % numero de subintervalos
N_cuad           = 250;
num_break_points = num_intervalos + 1; % numero de break points contando los extremos
kord             = 5; % orden de los bsplines, el grado de los polinomios es kord-1
regularity       = kord-2; % repeticion de los knots, simplemente degenerados por eso kord-2
xMax             = 5; % extremo superior del intervalo
xMin             = 0; % extremo inferior del intervalo
N_base           = num_intervalos + kord - 1;

%% armo la matriz de control points, si quiero toda la base c es la identidad
c = eye(N_base);

%% armamos los knots
[knots, zeta] = kntuniform (num_break_points, kord-1, regularity);

knots = (xMax - xMin)*knots + xMin;
zeta = (xMax - xMin)*zeta + xMin;

[x, w] = GaussLegendre_2(N_cuad);

x = 0.5*(xMax - xMin)*x + 0.5*(xMax + xMin);

%% calculo los bsplines en el vector x
bs = bspeval(kord-1, c, knots, x);

splines = [x'; bs]';
save('-ascii', 'splines.dat', 'splines')

%% para calcular las derivadas
[dc, dknots] = bspderiv(kord-1, c, knots);
[ddc, ddknots] = bspderiv(kord-2, dc, dknots);

dbs = bspeval(kord-2, dc, dknots, x);

ddbs = bspeval(kord-3, ddc, ddknots, x);

dsplines = [x'; dbs]';
save('-ascii', 'dsplines.dat', 'dsplines')

ddsplines = [x'; ddbs]';
save('-ascii', 'ddsplines.dat', 'ddsplines')

figure; hold on
for i = 1:N_base
  plot(splines(:,1), splines(:,i+1))
end
