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
Rc_i = 0/a0; % atomic units
Rc_f = 10/a0; % atomic units
N_dim_H = 8*N_base;
num_puntos_Rc = 20;

Rc_vec = linspace(Rc_i, Rc_f, num_puntos_Rc);

% Parametros de los materiales
% CdSe
m0_CdSe         = 0.13;
Eg_CdSe      = 1.714/eV; %1.751/eV; % Gap en atomic units
Ep_CdSe      = 24.27/eV; %17.4/eV; % atomic units
Delta_CdSe   = 0.24/eV; % 0.22/eV; % atomic units
gamma_CdSe   = [4.40, 1.6, 2.68]; %[5.51, 1.24, 2.14];
gamma_t_CdSe = gamma_CdSe;
% gamma_t_CdSe = [gamma_CdSe(1) - Ep_CdSe/(3*Eg_CdSe+Delta_CdSe),...
%                 gamma_CdSe(2) - 0.5*Ep_CdSe/(3*Eg_CdSe+Delta_CdSe),...
%                 gamma_CdSe(3) - 0.5*Ep_CdSe/(3*Eg_CdSe+Delta_CdSe)];

% ZnS (Al_0.5Ga_0.6As)
m0_ZnS       = 0.28;
Eg_ZnS       = 3.68/eV; %3.75/eV; % Gap en atomic units
Ep_ZnS       = 24.85/eV; %20.4/eV; % atomic units
Delta_ZnS    = 0.074/eV; % atomic units
gamma_ZnS    = [2.12, 0.51, 1.56]; %[1.77, 0.3, 0.62];
gamma_t_ZnS  = gamma_ZnS;
% gamma_t_ZnS  = [gamma_ZnS(1) - Ep_ZnS/(3*Eg_ZnS+Delta_ZnS),...
%                 gamma_ZnS(2) - 0.5*Ep_ZnS/(3*Eg_ZnS+Delta_ZnS),...
%                 gamma_ZnS(3) - 0.5*Ep_ZnS/(3*Eg_ZnS+Delta_ZnS)];

% offset
CBO = 0.9/eV; % atomic units (from J. App. Phys. 113 134304 (2013))
VBO = (Eg_ZnS - CBO - Eg_CdSe); %1.0990/eV; % atomic units

% Archivo de salida
name = sprintf('./resultados/E2.dat');

file = fopen(name, 'w');
fprintf(file, '# Num_interavaos = %i \n', num_intervalos);
fprintf(file, '# kord = %i \n', kord);
fprintf(file, '# N_base = %i , N_cuad = %i \n', N_base, N_cuad);
fprintf(file, '# N_dim hamiltoniano, N_dim_H = %i \n', N_dim_H);
fprintf(file, '# Parametros de CdSe \n');
fprintf(file, '# Eg = %f eV\n', eV*Eg_CdSe);
fprintf(file, '# Ep = %f eV\n', eV*Ep_CdSe);
fprintf(file, '# delta = %f eV\n', eV*Delta_CdSe);
fprintf(file, '# gamma = [%f, %f, %f] \n', gamma_CdSe);
fprintf(file, '# gamma_tilde = [%f, %f, %f] \n', gamma_t_CdSe);
fprintf(file, '# Parametros de ZnS \n');
fprintf(file, '# Eg = %f eV\n', eV*Eg_ZnS);
fprintf(file, '# Ep = %f eV\n', eV*Ep_ZnS);
fprintf(file, '# delta = %f eV\n', eV*Delta_ZnS);
fprintf(file, '# gamma = [%f, %f, %f] \n', gamma_ZnS);
fprintf(file, '# gamma_tilde = [%f, %f, %f] \n', gamma_t_ZnS);
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

for Rc = Rc_vec

  R1 = Rc + 0.8/a0;
  R2 = R1 + 3.5/a0;
  R3 = R2 + 1.0/a0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  me = m0_ZnS + (m0_CdSe - m0_ZnS)*heaviside(x-Rc).*heaviside(R1-x) ...
              + (m0_CdSe - m0_ZnS)*heaviside(x-R2).*heaviside(R3-x);

  Eg = Eg_ZnS + (Eg_CdSe - Eg_ZnS)*heaviside(x-Rc).*heaviside(R1-x) ...
              + (Eg_CdSe - Eg_ZnS)*heaviside(x-R2).*heaviside(R3-x);

  Ep = Ep_ZnS + (Ep_CdSe - Ep_ZnS)*heaviside(x-Rc).*heaviside(R1-x) ...
              + (Ep_CdSe - Ep_ZnS)*heaviside(x-R2).*heaviside(R3-x);

  Delta = Delta_ZnS + (Delta_CdSe - Delta_ZnS)*heaviside(x-Rc).*heaviside(R1-x) ...
                    + (Delta_CdSe - Delta_ZnS)*heaviside(x-R2).*heaviside(R3-x);

  gamma = gamma_ZnS + (gamma_CdSe - gamma_ZnS).*repmat(heaviside(x-Rc).*heaviside(R1-x),1,3) ...
                    + (gamma_CdSe - gamma_ZnS).*repmat(heaviside(x-R2).*heaviside(R3-x),1,3);

  gamma_t = gamma_t_ZnS + (gamma_t_CdSe - gamma_t_ZnS).*repmat(heaviside(x-Rc).*heaviside(R1-x),1,3) ...
                        + (gamma_t_CdSe - gamma_t_ZnS).*repmat(heaviside(x-R2).*heaviside(R3-x),1,3);

  Ev = 0.0 + VBO*heaviside(x-Rc).*heaviside(R1-x) ...
           + VBO*heaviside(x-R2).*heaviside(R3-x);

  Ec = Eg_ZnS + (VBO + Eg_CdSe - Eg_ZnS)*heaviside(x-Rc).*heaviside(R1-x) ...
              + (VBO + Eg_CdSe - Eg_ZnS)*heaviside(x-R2).*heaviside(R3-x);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculos las matrices que necesito para los operadores
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Primero calculo la matriz de solapamiento
  S = bs*(repmat(w.*x.^2, 1, N_base).*bs');
  S = triu(S) + triu(S, 1)';

  Ev_matrix = bs*(repmat(w.*x.^2.*Ev, 1, N_base).*bs');
  Ec_matrix = bs*(repmat(w.*x.^2.*Ec, 1, N_base).*bs');

  A = Ec_matrix + hbar^2/(2*m0)*dbs*(repmat(w.*x.^2.*me, 1, N_base).*dbs');

  P = -Ev_matrix + hbar^2/(2*m0)*dbs*(repmat(w.*x.^2.*gamma_t(:,1).*me, 1, N_base).*dbs');

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

  [vec, auval] = eigs(HT, ST, size(HT, 1)-2);
  [auval, I] = sort(real(eV*diag(auval)));

  auvec = vec(:,I);

  fprintf(file, '%17.9e  ', Rc*a0);
  for ind_auval = 1:size(auval)
    fprintf(file, '%17.9e  ', auval(ind_auval));
  end
  fprintf(file, '\n');

end
