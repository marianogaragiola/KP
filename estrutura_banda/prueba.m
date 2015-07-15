clear all;
close all;

% Parametros independientes de los materiales.
m0 = 0.510998928E6; % masa del electron en eV/C^2
h = 6.58211899E-16; % hbarra en eV/S


gamma_L_1 = [19.67, 8.37, 9.29];
E_c = 0.418; % en eV
E_v = 0.0; % en eV
Ep = 22.2; % en eV
Delta = 0.38; % en eV

% Parametros dependientes.
Eg = E_c-E_v; % en eV
P0 = h^2*Ep/(2.0*m0);

% vector gamma
gamma(1) = gamma_L(1) - Ep/(3.0*Eg+Delta);
gamma(2:3) = gamma_L(2:3) - 0.5*Ep/(3.0*Eg+Delta);

% vecto k
k = [0, 0, 1];

% Calculo los elementos de matriz
A = function_A(h, m0, E_c, E_v, P0, gamma, k);
P = function_P(h, m0, E_c, E_v, P0, gamma, k);
Q = function_Q(h, m0, E_c, E_v, P0, gamma, k);
R = function_R(h, m0, E_c, E_v, P0, gamma, k);
S = function_S(h, m0, E_c, E_v, P0, gamma, k);
U = function_S(h, m0, E_c, E_v, P0, gamma, k);
V = function_V(h, m0, E_c, E_v, P0, gamma, k);

% Armo el hamiltoniano completo como una matriz 8x8
H = [A, 0, conj(V), 0, sqrt(3.)*V, -sqrt(2.)*U, -U, sqrt(2.)*conj(V);
     0, A, -sqrt(2.)*U, -sqrt(3.)*conj(V), 0, -V, sqrt(2.)*V, U;
     V, -sqrt(2.)*U, -P+Q, -conj(S), R, 0, sqrt(3./2.)*S, -sqrt(2.)*Q;
     0, -sqrt(3.)*V, -S, -P-Q, 0, R, -sqrt(2.)*R, 1./sqrt(2.)*S;
     sqrt(3.)*conj(V), 0, conj(R), 0, -P-Q, conj(S), 1./sqrt(2.)*conj(S), sqrt(2)*conj(R);
     -sqrt(2.)*U, -conj(V), 0, conj(R), S, -P+Q, sqrt(2.)*Q, sqrt(3./2.)*conj(S);
     -U, sqrt(2.)*conj(V), sqrt(3./2.)*conj(S), -sqrt(2.)*R, 1./sqrt(2.)*S, sqrt(2.)*Q, -P-Delta, 0;
     sqrt(2.)*V, U, -sqrt(2.)*Q, 1./sqrt(2.)*conj(S), sqrt(2.)*R, sqrt(3./2.)*S, 0, -P-Delta];
