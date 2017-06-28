clear all
close all

format long
pkg load nurbs

% Parametros para los B-splines
num_intervalos   = 50; % numero de subintervalos
N_cuad           = 100;
num_break_points = num_intervalos + 1; % numero de break points contando los extremos
kord             = 7; % orden de los bsplines, el grado de los polinomios es kord-1
regularity       = kord-2; % repeticion de los knots, simplemente degenerados por eso kord-2
xMax             = 20; % extremo superior del intervalo
xMin             = 0; % extremo inferior del intervalo
N_splines        = num_intervalos + kord - 1;
N_base           = N_splines - 1;
num_knots        = N_splines + kord;

% Parametros fisicos.
m0 = 0.510998928E6; % masa del electron en eV/C^2
m0 = 0.510998928E6/1E18; % masa del electron en eV*s^2/nm^2
hbar = 6.58211899E-16; % hbarra en eV*S
m_angular = 0; % componente en el eje z del momento angular.
k_z = 0; % Componente del vector k en el eje z

% Parametros del material (ahora uso de GaAs)
Eg = 1.519; % [eV]
Ep = 23.81; % [eV]
delta = 0.341; % [eV]
E1 = Eg; % Bottom of the conduction band
E2 = 0.0; % Top of the valence band
P = sqrt(0.5*Ep*hbar^2/m0);
gamma = [7.05 2.35 3];% [1.8251, -0.2625, 0.3875]; % Son los \tilde{gamma} del paper de Peeters

%% armo la matriz de control points, si quiero toda la base c es la identidad
c = eye(N_splines);
c(N_splines, N_splines) = 0;

%% armamos los knots
[knots, zeta] = kntuniform (num_break_points, kord-1, regularity);

knots = (xMax - xMin)*knots + xMin;
zeta = (xMax - xMin)*zeta + xMin;

[x, w] = GaussLegendre_2(N_cuad);
% x = 0.5*(xMax - xMin)*x + 0.5*(xMax + xMin);

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

% Calculos las matrices.
% Primero calculo la matriz de solapamiento
S = bs*(repmat(w.*x, 1, N_base).*bs');

% Caculo la matriz de derivadas
T = dbs*(repmat(w.*x, 1, N_base).*dbs');

% Calculo la matriz de B_{n1} dB_{n2}/dr
D = bs*(repmat(w.*x, 1, N_base).*dbs');

% Matriz del operador 1/rho^2
rho_2 = bs*(repmat(w./x, 1, N_base).*bs');

% Matriz del operador 1/rho
rho_1 = bs*(repmat(w, 1, N_base).*bs');

% Calculo el bloque E_{7-}^H
E7MH = E1*S + 0.5*hbar^2/m0*( T + m_angular^2*rho_2 + k_z^2*S);

% Calculo el bloque E_{7+}^H
E7PH = (E2 - delta)*S - 0.5*hbar^2/m0*gamma(1)*(T + m_angular^2*rho_2 + k_z^2*S);

% Calculo el bloque E_{8+}^H
E8PH = E2*S - 0.5*hbar^2/m0*((gamma(1)+gamma(2))*(T + m_angular^2*rho_2) + ...
                             (gamma(1)-2*gamma(2))*k_z^2*S);

% Calculo el bloque E_{8-}^H
E8PL = E2*S - 0.5*hbar^2/m0*((gamma(1)-gamma(2))*(T + m_angular^2*rho_2) + ...
                             (gamma(1)+2*gamma(2))*k_z^2*S);

% Calculo el bloque de P^+
% PP = -i*P*kronDel(m_angular2-m_angular1+1)*(D - m_angular2*rho_1);

% Calculo el bloque de P^z
Pz = P*k_z*S;

% Calculo el bloque A
A = 0.5*hbar^2*gamma(2)/m0*(-T - m_angular^2*rho_2 - 2.0*k_z^2*S);

% Faltan los bloques de B y C pero ellos cambian el momento angular
% Tengo que averiguar como se hace eso.


% Armo la matriz del hamiltoniano
H = speye(8*N_base);
% Armo los bloques diagonales del hamiltoniano
H(1:N_base,1:N_base)                       = E7MH(:,:);
H(N_base+1:2*N_base,N_base+1:2*N_base)     = E7MH(:,:);
H(2*N_base+1:3*N_base,2*N_base+1:3*N_base) = E8PH(:,:);
H(3*N_base+1:4*N_base,3*N_base+1:4*N_base) = E8PL(:,:);
H(4*N_base+1:5*N_base,4*N_base+1:5*N_base) = E8PL(:,:);
H(5*N_base+1:6*N_base,5*N_base+1:6*N_base) = E8PH(:,:);
H(6*N_base+1:7*N_base,6*N_base+1:7*N_base) = E7PH(:,:);
H(7*N_base+1:8*N_base,7*N_base+1:8*N_base) = E7PH(:,:);

H(3*N_base+1:4*N_base,6*N_base+1:7*N_base) = -sqrt(2.0)*A(:,:);
H(6*N_base+1:7*N_base,3*N_base+1:4*N_base) = -sqrt(2.0)*A(:,:);
H(4*N_base+1:5*N_base,7*N_base+1:8*N_base) = sqrt(2.0)*A(:,:);
H(7*N_base+1:8*N_base,4*N_base+1:5*N_base) = sqrt(2.0)*A(:,:);

H(1:N_base,3*N_base+1:4*N_base)            = sqrt(2./3.)*Pz(:,:);
H(3*N_base+1:4*N_base,1:N_base)            = sqrt(2./3.)*Pz(:,:);
H(1:N_base,6*N_base+1:7*N_base)            = sqrt(1./3.)*Pz(:,:);
H(6*N_base+1:7*N_base,1:N_base)            = sqrt(1./3.)*Pz(:,:);
H(N_base+1:2*N_base,4*N_base+1:5*N_base)   = sqrt(2./3.)*Pz(:,:);
H(4*N_base+1:5*N_base,N_base+1:2*N_base)   = sqrt(2./3.)*Pz(:,:);
H(N_base+1:2*N_base,7*N_base+1:8*N_base)   = -sqrt(1./3.)*Pz(:,:);
H(7*N_base+1:8*N_base,N_base+1:2*N_base)   = -sqrt(1./3.)*Pz(:,:);

ST = speye(8*N_base);
for i = 1:8
  ST((i-1)*N_base+1:i*N_base,(i-1)*N_base+1:i*N_base) = S(:,:);
end

auval = eig(H, ST);
auval = sort(auval);

plot(auval,'-o')
