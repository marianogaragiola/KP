clear all
close all

format long
addpath(genpath('~/octave/nurbs-1.3.13/'));

a0 = 0.0529177210; eV = 27.21138564; c_light = 137.035999074492;
alpha = 658.4092645439;

% Parametros para los B-splines
num_intervalos   = 100; % numero de subintervalos
N_cuad           = 100;
num_break_points = num_intervalos + 1; % numero de break points contando los extremos
kord             = 7; % orden de los bsplines, el grado de los polinomios es kord-1
regularity       = kord-2; % repeticion de los knots, simplemente degenerados por eso kord-2
RMax             = 20/a0; % extremo superior del intervalo
RMin             = 0; % extremo inferior del intervalo
N_splines        = num_intervalos + kord - 1;
N_base           = N_splines - 1;
num_knots        = N_splines + kord;

% Parametros fisicos.
m0 = 1.0; % electron mass in atomic units
hbar = 1.0; % Planck constant in atomic units.
m_angular = 0; % Momento angular total
mu_B = 5.7883818066e-5; % Borh magneton in eV/T
R1_i = 0/a0; % atomic units
R1_f = 14/a0; % atomic units
R2 = 15/a0; % atomic units
N_dim_H = 8*(2*m_angular+1)*N_base;
num_puntos_R1 = 50;

R1_vec = linspace(R1_i, R1_f, num_puntos_R1);

% Parametros de los materiales
% GaAs
Eg_GaAs      = 1.519/eV; % Gap en atomic units
Ep_GaAs      = 23.81/eV; % atomic units
Delta_GaAs   = 0.341/eV; % atomic units
gamma_GaAs   = [7.05, 2.35, 3.0];
gamma_t_GaAs = [1.8251, -0.2625, 0.3875];

% AlGaAs (Al_0.5Ga_0.6As)
Eg_AlGaAs       = 2.1634/eV; % Gap en atomic units
Ep_AlGaAs       = 22.726/eV; % atomic units
Delta_AlGaAs    = 0.3166/eV; % atomic units
gamma_AlGaAs    = [5.1217, 1.4347, 2.0971];
gamma_t_AlGaAs  = [1.6201, -0.3161, 0.3463];

% offset
CBO = 0.348/eV; % atomic units
VBO = 0.228/eV; % atomic units

% Archivo de salida
name = sprintf('./resultados/E.dat');

file = fopen(name, 'w');
fprintf(file, '# Num_interavaos = %i \n', num_intervalos);
fprintf(file, '# kord = %i \n', kord);
fprintf(file, '# N_base = %i , N_cuad = %i \n', N_base, N_cuad);
fprintf(file, '# N_dim hamiltoniano, N_dim_H = %i \n', N_dim_H);
fprintf(file, '# Parametros de GaAs \n');
fprintf(file, '# Eg = %f eV\n', eV*Eg_GaAs);
fprintf(file, '# Ep = %f eV\n', eV*Ep_GaAs);
fprintf(file, '# delta = %f eV\n', eV*Delta_GaAs);
fprintf(file, '# gamma = [%f, %f, %f] \n', gamma_GaAs);
fprintf(file, '# gamma_tilde = [%f, %f, %f] \n', gamma_t_GaAs);
fprintf(file, '# Parametros de AlGaAs \n');
fprintf(file, '# Eg = %f eV\n', eV*Eg_AlGaAs);
fprintf(file, '# Ep = %f eV\n', eV*Ep_AlGaAs);
fprintf(file, '# delta = %f eV\n', eV*Delta_AlGaAs);
fprintf(file, '# gamma = [%f, %f, %f] \n', gamma_AlGaAs);
fprintf(file, '# gamma_tilde = [%f, %f, %f] \n', gamma_t_AlGaAs);
fprintf(file, '# CBO = %f, VBO = %f\n', eV*CBO, eV*VBO);
fprintf(file, '# autovalores del problema \n');

%% armo la matriz de control points, si quiero toda la base c es la identidad
c = eye(N_splines);
c(N_splines, N_splines) = 0;

%% armamos los knots
[knots, zeta] = kntuniform (num_break_points, kord-1, regularity);

knots = (RMax - RMin)*knots + RMin;
zeta = (RMax - RMin)*zeta + RMin;

[x, w] = GaussLegendre_2(N_cuad);

aux1 = []; aux2 = [];
for i = kord+1:num_knots-kord+1
  aux1 = [aux1; 0.5*(knots(i) - knots(i-1))*x + 0.5*(knots(i) + knots(i-1))];
  aux2 = [aux2; 0.5*(knots(i) - knots(i-1))*w];
end

x = aux1;
w = aux2;

%% calculo los bsplines en el vector x
bs = bspeval(kord-1, c, knots, x);
bs = sparse(bs(1:N_base,:));

%% para calcular las derivadas
[dc, dknots] = bspderiv(kord-1, c, knots);

% calculo las deridabas de los bsplines
dbs = bspeval(kord-2, dc, dknots, x);
dbs = sparse(dbs(1:N_base,:));

for R1 = R1_vec
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Eg = Eg_AlGaAs + (Eg_GaAs - Eg_AlGaAs)*heaviside(x-R1).*heaviside(R2-x);
    Ep = Ep_AlGaAs + (Ep_GaAs - Ep_AlGaAs)*heaviside(x-R1).*heaviside(R2-x);
    Delta = Delta_AlGaAs + (Delta_GaAs - Delta_AlGaAs)*heaviside(x-R1).*heaviside(R2-x);
    gamma = gamma_AlGaAs + (gamma_GaAs - gamma_AlGaAs).*repmat(heaviside(x-R1).*heaviside(R2-x),1,3);
    gamma_t = gamma_t_AlGaAs + (gamma_t_GaAs - gamma_t_AlGaAs).*repmat(heaviside(x-R1).*heaviside(R2-x),1,3);

    Ev = 0.0 + VBO*heaviside(x-R1).*heaviside(R2-x);
    Ec = Eg_AlGaAs + (VBO + Eg_GaAs - Eg_AlGaAs)*heaviside(x-R1).*heaviside(R2-x);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculos las matrices que necesito para los operadores
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Primero calculo la matriz de solapamiento
    S = bs*(repmat(w.*x.^2, 1, N_base).*bs');
    S = triu(S) + triu(S, 1)';

    Ev_matrix = bs*(repmat(w.*x.^2.*Ev, 1, N_base).*bs');
    Ec_matrix = bs*(repmat(w.*x.^2.*Ec, 1, N_base).*bs');

    A = Ec_matrix + hbar^2/(2*m0)*dbs*(repmat(w.*x.^2, 1, N_base).*dbs');

    P = -Ev_matrix + hbar^2/(2*m0)*dbs*(repmat(w.*x.^2.*gamma_t(:,1), 1, N_base).*dbs');

    delta_matrix = bs*(repmat(w.*x.^2.*Delta, 1, N_base).*bs');

    HT = speye(8*N_base);
    ST = speye(8*N_base);

    HT(0*N_base+1:1*N_base,0*N_base+1:1*N_base) = A; % Bloque diagonal
    HT(1*N_base+1:2*N_base,1*N_base+1:2*N_base) = A;
    HT(2*N_base+1:3*N_base,2*N_base+1:3*N_base) = -P;
    HT(3*N_base+1:4*N_base,3*N_base+1:4*N_base) = -P;
    HT(4*N_base+1:5*N_base,4*N_base+1:5*N_base) = -P;
    HT(5*N_base+1:6*N_base,5*N_base+1:6*N_base) = -P;
    HT(6*N_base+1:7*N_base,6*N_base+1:7*N_base) = -P-delta_matrix;
    HT(7*N_base+1:8*N_base,7*N_base+1:8*N_base) = -P-delta_matrix;

    ST(0*N_base+1:1*N_base,0*N_base+1:1*N_base) = S; % Bloque diagonal
    ST(1*N_base+1:2*N_base,1*N_base+1:2*N_base) = S;
    ST(2*N_base+1:3*N_base,2*N_base+1:3*N_base) = S;
    ST(3*N_base+1:4*N_base,3*N_base+1:4*N_base) = S;
    ST(4*N_base+1:5*N_base,4*N_base+1:5*N_base) = S;
    ST(5*N_base+1:6*N_base,5*N_base+1:6*N_base) = S;
    ST(6*N_base+1:7*N_base,6*N_base+1:7*N_base) = S;
    ST(7*N_base+1:8*N_base,7*N_base+1:8*N_base) = S;

    auval = eigs(HT, ST, size(HT, 1)-2);
    auval = eV*sort(real(auval));

    fprintf(file, '%17.9e  ', R1*a0);
    for ind_auval = 1:size(auval)
      fprintf(file, '%17.9e  ', auval(ind_auval));
    end
    fprintf(file, '\n');

end
