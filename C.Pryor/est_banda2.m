clear all;
close all;

% Parametros independientes de los materiales.
m0 = 0.510998928E6; % masa del electron en eV/C^2
h = 6.58211899E-16; % hbarra en eV/S

% Radio del core del quantum dot
Rc = 10; % en nm
radios = linspace(1, 20, 100);

% Parametros del material del Shell (GaAs en este caso)
gamma_L_2 = [6.85, 2.1, 2.9];
Eg_2 = 1.519;
E_v_2 = 0.0; % en eV
E_c_2 = Eg_2 + E_v_2; % en eV
Ep_2 = 25.7; % en eV
Delta_2 = 0.33; % en eV

% Parametros del material del Core (InAs en este caso)
gamma_L_1 = [19.67, 8.37, 9.29];
Eg_1 = 0.418; % en eV
V0 = 0.085; % en eV offset entre las bandas de valencias del Core y el Shell
E_v_1 = E_v_2 + V0; % en eV
E_c_1 = Eg_1 + E_v_1; % en eV
Ep_1 = 22.2; % en eV
Delta_1 = 0.38; % en eV

salida = [];

for r = radios

  % Parametros dependientes.
  Eg = Eg_2 + (Eg_1-Eg_2)*(1.0-heaviside(r-Rc));
  Ep = Ep_2 + (Ep_1-Ep_2)*(1.0-heaviside(r-Rc));
  E_c = E_c_2 + (E_c_1-E_c_2)*(1.0-heaviside(r-Rc));
  E_v = E_v_2 + (E_v_1-E_v_2)*(1.0-heaviside(r-Rc));; % en eV
  Delta = Delta_2 + (Delta_1-Delta_2)*(1.0-heaviside(r-Rc));
  gamma_L = gamma_L_2 + (gamma_L_1-gamma_L_2)*(1.0-heaviside(r-Rc));

  P0 = h^2*Ep/(2.0*m0);

  % vector gamma
  gamma(1) = gamma_L(1) - Ep/(3.0*Eg+Delta);
  gamma(2:3) = gamma_L(2:3) - 0.5*Ep/(3.0*Eg+Delta);

  % vecto k
  k = [1, 1, 1];
  k = k./norm(k);

  % Calculo los elementos de matriz
  A = function_A(h, m0, E_c, E_v, P0, gamma, k);
  P = function_P(h, m0, E_c, E_v, P0, gamma, k);
  Q = function_Q(h, m0, E_c, E_v, P0, gamma, k);
  R = function_R(h, m0, E_c, E_v, P0, gamma, k);
  S = function_S(h, m0, E_c, E_v, P0, gamma, k);
  U = function_U(h, m0, E_c, E_v, P0, gamma, k);
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

  e = eig(H);
  e = e';

  salida = [salida;r e];

end

save('./resultados/k111.dat', '-ascii', 'salida');
