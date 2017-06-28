% Codigo que calcula el hamiltoniano kp segun
% J. Phys.: Condens. Matter 26 (2014) 095501
% tiene en cuenta el momento angular guarda
% los autovalores en funcion de k_z.
%
% Usa B-splines como funciones base

clear all
close all

format long
addpath(genpath('~/octave/nurbs-1.3.13/'));

a0 = 0.0529177210; eV = 27.21138564; c_light = 137.035999074492;
alpha = 658.4092645439;

% Parametros para los B-splines
num_intervalos   = 50; % numero de subintervalos
N_cuad           = 100;
num_break_points = num_intervalos + 1; % numero de break points contando los extremos
kord             = 7; % orden de los bsplines, el grado de los polinomios es kord-1
regularity       = kord-2; % repeticion de los knots, simplemente degenerados por eso kord-2
xMax             = 10/a0; % extremo superior del intervalo
xMin             = 0; % extremo inferior del intervalo
N_splines        = num_intervalos + kord - 1;
N_base           = N_splines - 1;
num_knots        = N_splines + kord;

% Parametros fisicos.
m0 = 1.0; % electron mass in atomic units
hbar = 1.0; % Planck constant in atomic units.
m_angular = 4; % Momento angular total
mu_B = 5.7883818066e-5; % Borh magneton in eV/T
kz_i = 0.*a0; % Componente del vector k en el eje z
kz_f = 0.8*a0; %
num_puntos_kz = 40;

N_dim_H = 8*(2*m_angular+1)*N_base;

% Parametros del material (ahora uso de GaAs)
Eg = 1.519/eV; % atomic units
Ep = 23.81/eV; % atomic units
EPX = 15.79/eV; % atomic units
EGC = 0;
delta = 0.341/eV; % atomic units
delta_C = 0.180/eV; % atomic units
mu_B = mu_B/eV; % paso atomic units
E1 = Eg; % Bottom of the conduction band
E2 = 0.0; % Top of the valence band
P = sqrt(0.5*Ep*hbar^2/m0);

% Parametros de Peeters
gamma    = [7.05, 2.35, 3];%
gamma_t  = [1.8251, -0.2625, 0.3875]; % Son los \tilde{gamma} del paper de Peeters
gamma_d  = [5.8967, 2.1516, 2.7194];
gamma_td = [1.6297, -0.2214, 0.3465];

gamma = gamma_t;
gamma_d = gamma_td;

kappa = [1, -1, 3, 1, -1, -3, 2, -2];

g_CB = 2 - 2/3*Eg*delta/(Eg*(Eg + delta)); % Lande factor for the conduction band
g_VB = gamma(3) + 2/3*gamma(2) - 1/3*gamma(1) - 2/3;  % Lande factor for the valence band

kz_vec = linspace(kz_i, kz_f, num_puntos_kz);

% Archivo de salida
name = sprintf('./resultados/E_vs_kz.dat');

file = fopen(name, 'w');
fprintf(file, '# Num_interavaos = %i \n', num_intervalos);
fprintf(file, '# kord = %i \n', kord);
fprintf(file, '# N_base = %i , N_cuad = %i \n', N_base, N_cuad);
fprintf(file, '# N_dim hamiltoniano, N_dim_H = %i \n', N_dim_H);
fprintf(file, '# Parametros del material \n');
fprintf(file, '# Eg = %f eV\n', eV*Eg);
fprintf(file, '# Ep = %f eV\n', eV*Ep);
fprintf(file, '# Epx = %f eV\n', eV*EPX);
fprintf(file, '# EGC = %f eV\n', eV*EGC);
fprintf(file, '# delta = %f eV\n', eV*delta);
fprintf(file, '# Delta_c = %f eV\n' , eV*delta_C);
fprintf(file, '# gamma = [%f, %f, %f] \n', gamma);
fprintf(file, '# gamma_tilde = [%f, %f, %f] \n', gamma_t);
fprintf(file, '# gamma_delta = [%f, %f, %f] \n', gamma_d);
fprintf(file, '# gamma_td = [%f, %f, %f] \n', gamma_td);
fprintf(file, '# kz inicial = %f T\n', kz_i);
fprintf(file, '# kz final = %f T\n', kz_f);
fprintf(file, '# autovalores del problema \n');

%% armo la matriz de control points, si quiero toda la base c es la identidad
c = eye(N_splines);
c(N_splines, N_splines) = 0;

%% armamos los knots
[knots, zeta] = kntuniform (num_break_points, kord-1, regularity);

knots = (xMax - xMin)*knots + xMin;
zeta = (xMax - xMin)*zeta + xMin;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculos las matrices que necesito para los operadores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primero calculo la matriz de solapamiento
S = bs*(repmat(w.*x, 1, N_base).*bs');
S = triu(S) + triu(S, 1)';

% Caculo la matriz de derivadas
T = dbs*(repmat(w.*x, 1, N_base).*dbs');
T = triu(T) + triu(T, 1)';

% Calculo la matriz de rho B_{n1} dB_{n2}/dr
D = bs*(repmat(w.*x, 1, N_base).*dbs');
% D = triu(D) + triu(D, 1)';

% Calculo la matriz de B_{n1} dB_{n2}/dr
BdB = bs*(repmat(w, 1, N_base).*dbs');
% S = triu(S) + triu(S, 1)';

% Calculo la matriz de rho B_{n1} dB_{n2}/dr
D_rho = bs*(repmat(w.*x.^2,1, N_base).*dbs');
% S = triu(S) + triu(S, 1)';

% Matriz del operador rho
rho1 = bs*(repmat(w.*x.^2, 1, N_base).*bs');
rho1 = triu(rho1) + triu(rho1, 1)';

% Matriz del operador rho^2
rho2 = bs*(repmat(w.*x.^3, 1, N_base).*bs');
rho2 = triu(rho2) + triu(rho2, 1)';

% Matriz del operador 1/rho^2
rho_2 = bs*(repmat(w./x, 1, N_base).*bs');
rho_2 = triu(rho_2) + triu(rho_2, 1)';

% Matriz del operador 1/rho
rho_1 = bs*(repmat(w, 1, N_base).*bs');
rho_1 = triu(rho_1) + triu(rho_1, 1)';

%%% LOOP EN EL CAMPO MAGNETICO
for ind_kz = 1:num_puntos_kz

  k_z = kz_vec(ind_kz);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Armo el hamiltoniano completo
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  HT = speye(8*(2*m_angular+1)*N_base);
  ST = speye(8*(2*m_angular+1)*N_base);

  for m1 = -m_angular:m_angular
    for m2 = -m_angular:m_angular

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Calculos las matrices de los operadores
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Calculo la matrix del operador k_{-}k_{+} = k_{+}k_{-}
      kpkm = (T + m2^2*rho_2);

      % Calculo la matriz del operador k_{-}k_{-}
      kmkm = kronDel(m2-m1-2,0)*(T + (2 - 2*m2)*BdB - m2*(m2-2)*rho_2);
      % kmkm = kronDel(m2-m1-2,0)*(T + (2 + 2*m2)*BdB - m2*(m2+2)*rho_2);

      % Calculo la matriz del operador k_{+}k_{+}
      kpkp = kronDel(m2-m1+2,0)*(T + (2 + 2*m2)*BdB - m2*(m2+2)*rho_2);
      % kpkp = kronDel(m2-m1+2,0)*(T + (2 - 2*m2)*BdB - m2*(m2-2)*rho_2);

      % Calculo el bloque E_{7-}^H
      E7MH = kronDel(m1,m2)*(E1*S + 0.5*hbar^2/m0*(kpkm + k_z^2*S));

      % Calculo el bloque E_{7+}^H
      % E7PH = kronDel(m1,m2)*((E2 - delta)*S - 0.5*hbar^2/m0*gamma(1)*(kpkm + k_z^2*S));
      E7PH = kronDel(m1,m2)*((E2 + delta)*S - 0.5*hbar^2/m0*gamma(1)*(kpkm + k_z^2*S));

      % Calculo el bloque E_{8+}^H
      E8PH = kronDel(m1,m2)*(E2*S - 0.5*hbar^2/m0*((gamma(1) + gamma(2))*kpkm + ...
                                                   (gamma(1) - 2*gamma(2))*k_z^2*S));

      % Calculo el bloque E_{8-}^H
      E8PL = kronDel(m1,m2)*(E2*S - 0.5*hbar^2/m0*((gamma(1) - gamma(2))*kpkm + ...
                                                   (gamma(1) + 2*gamma(2))*k_z^2*S));

      % Calculo el bloque A
      A = kronDel(m1,m2)*(0.5*hbar^2*gamma(2)/m0*(-kpkm - 2.0*k_z^2*S));
      % Calculo el bloque A_{\Delta}
      A_d = kronDel(m1,m2)*(0.5*hbar^2*gamma_d(2)/m0*(-kpkm - 2.0*k_z^2*S));

      % Calculo el bloque de P^+
      PP = -full((D - m2*rho_1));
      PP = P*kronDel(m2-m1+1,0)*sparse(complex(0, PP));

      % Calculo el bloque de P^-
      PM = -full((D + m2*rho_1));
      PM = P*kronDel(m2-m1-1,0)*sparse(complex(0, PM));

      % Calculo el bloque de P^z
      Pz = kronDel(m1,m2)*P*k_z*S;

      % Calculo el bloque B
      B = full(D + m2*rho_1);
      B = sqrt(3)*hbar^2*gamma(3)/m0*k_z*kronDel(m2-m1-1,0)*...
          sparse(complex(B, 0));

      % Calculo el bloque B_{\Delta}
      B_d = gamma_d(3)/gamma(3)*B;

      % Calculo el bloque C
      C = 0.25*sqrt(3)*hbar^2/m0*((gamma(2)+gamma(3))*kmkm + (gamma(2)+gamma(3))*kpkp);
      % C = 0.25*sqrt(3)*hbar^2/m0*((gamma(2)+gamma(3))*kmkm + (gamma(2)-gamma(3))*kpkp);

      % Calculo el bloque C_{\Delta}
      C_d = 0.25*sqrt(3)*hbar^2/m0*((gamma_d(2)+gamma_d(3))*kmkm + (gamma_d(2)+gamma_d(3))*kpkp);

      % % Armo la matriz del hamiltoniano
      H8x8 = speye(8*N_base);
      % Primer fila
      H8x8(0*N_base+1:1*N_base,0*N_base+1:1*N_base) = E7MH(:,:); % Bloque diagonal
      H8x8(0*N_base+1:1*N_base,1*N_base+1:2*N_base) = 0;
      H8x8(0*N_base+1:1*N_base,2*N_base+1:3*N_base) = -sqrt(1/2)*PP;
      H8x8(0*N_base+1:1*N_base,3*N_base+1:4*N_base) = sqrt(2/3)*Pz;
      H8x8(0*N_base+1:1*N_base,4*N_base+1:5*N_base) = sqrt(1/6)*PM;
      H8x8(0*N_base+1:1*N_base,5*N_base+1:6*N_base) = 0;
      H8x8(0*N_base+1:1*N_base,6*N_base+1:7*N_base) = sqrt(1/3)*Pz;
      H8x8(0*N_base+1:1*N_base,7*N_base+1:8*N_base) = sqrt(1/3)*PM;
      % Segunda fila
      H8x8(1*N_base+1:2*N_base,0*N_base+1:1*N_base) = 0;
      H8x8(1*N_base+1:2*N_base,1*N_base+1:2*N_base) = E7MH(:,:); % Bloque diagonal
      H8x8(1*N_base+1:2*N_base,2*N_base+1:3*N_base) = 0;
      H8x8(1*N_base+1:2*N_base,3*N_base+1:4*N_base) = -sqrt(1/6)*PP;
      H8x8(1*N_base+1:2*N_base,4*N_base+1:5*N_base) = sqrt(2/3)*Pz;
      H8x8(1*N_base+1:2*N_base,5*N_base+1:6*N_base) = sqrt(1/2)*PM;
      H8x8(1*N_base+1:2*N_base,6*N_base+1:7*N_base) = sqrt(1/3)*PP;
      H8x8(1*N_base+1:2*N_base,7*N_base+1:8*N_base) = -sqrt(1/3)*Pz;
      % Tercer fila
      H8x8(2*N_base+1:3*N_base,0*N_base+1:1*N_base) = -sqrt(1/2)*PM;
      H8x8(2*N_base+1:3*N_base,1*N_base+1:2*N_base) = 0;
      H8x8(2*N_base+1:3*N_base,2*N_base+1:3*N_base) = E8PH(:,:); % Bloque diagonal
      H8x8(2*N_base+1:3*N_base,3*N_base+1:4*N_base) = B;
      H8x8(2*N_base+1:3*N_base,4*N_base+1:5*N_base) = C;
      H8x8(2*N_base+1:3*N_base,5*N_base+1:6*N_base) = 0;
      H8x8(2*N_base+1:3*N_base,6*N_base+1:7*N_base) = sqrt(1/2)*B_d;
      H8x8(2*N_base+1:3*N_base,7*N_base+1:8*N_base) = sqrt(2)*C_d;
      % Cuarta fila
      H8x8(3*N_base+1:4*N_base,0*N_base+1:1*N_base) = sqrt(2/3)*Pz;
      H8x8(3*N_base+1:4*N_base,1*N_base+1:2*N_base) = -sqrt(1/6)*PM;
      H8x8(3*N_base+1:4*N_base,2*N_base+1:3*N_base) = B';
      H8x8(3*N_base+1:4*N_base,3*N_base+1:4*N_base) = E8PL(:,:); % Bloque diagonal
      H8x8(3*N_base+1:4*N_base,4*N_base+1:5*N_base) = 0;
      H8x8(3*N_base+1:4*N_base,5*N_base+1:6*N_base) = C;
      H8x8(3*N_base+1:4*N_base,6*N_base+1:7*N_base) = -sqrt(2)*A_d;
      H8x8(3*N_base+1:4*N_base,7*N_base+1:8*N_base) = -sqrt(3/2)*B_d;
      % Quinta fila
      H8x8(4*N_base+1:5*N_base,0*N_base+1:1*N_base) = sqrt(1/6)*PP;
      H8x8(4*N_base+1:5*N_base,1*N_base+1:2*N_base) = sqrt(2/3)*Pz;
      H8x8(4*N_base+1:5*N_base,2*N_base+1:3*N_base) = C';
      H8x8(4*N_base+1:5*N_base,3*N_base+1:4*N_base) = 0;
      H8x8(4*N_base+1:5*N_base,4*N_base+1:5*N_base) = E8PL(:,:); % Bloque diagonal
      H8x8(4*N_base+1:5*N_base,5*N_base+1:6*N_base) = -B;
      H8x8(4*N_base+1:5*N_base,6*N_base+1:7*N_base) = -sqrt(3/2)*B_d';
      H8x8(4*N_base+1:5*N_base,7*N_base+1:8*N_base) = sqrt(2)*A_d;
      % Sexta fila
      H8x8(5*N_base+1:6*N_base,0*N_base+1:1*N_base) = 0;
      H8x8(5*N_base+1:6*N_base,1*N_base+1:2*N_base) = sqrt(1/2)*PP;
      H8x8(5*N_base+1:6*N_base,2*N_base+1:3*N_base) = 0;
      H8x8(5*N_base+1:6*N_base,3*N_base+1:4*N_base) = C';
      H8x8(5*N_base+1:6*N_base,4*N_base+1:5*N_base) = -B';
      H8x8(5*N_base+1:6*N_base,5*N_base+1:6*N_base) = E8PH(:,:); % Bloque diagonal
      H8x8(5*N_base+1:6*N_base,6*N_base+1:7*N_base) = -sqrt(2)*C_d';
      H8x8(5*N_base+1:6*N_base,7*N_base+1:8*N_base) = sqrt(1/2)*B_d';
      % Septima fila
      H8x8(6*N_base+1:7*N_base,0*N_base+1:1*N_base) = sqrt(1/3)*Pz;
      H8x8(6*N_base+1:7*N_base,1*N_base+1:2*N_base) = sqrt(1/3)*PM;
      H8x8(6*N_base+1:7*N_base,2*N_base+1:3*N_base) = sqrt(1/2)*B_d';
      H8x8(6*N_base+1:7*N_base,3*N_base+1:4*N_base) = -sqrt(2)*A_d';
      H8x8(6*N_base+1:7*N_base,4*N_base+1:5*N_base) = -sqrt(3/2)*B_d;
      H8x8(6*N_base+1:7*N_base,5*N_base+1:6*N_base) = -sqrt(2)*C_d;
      H8x8(6*N_base+1:7*N_base,6*N_base+1:7*N_base) = E7PH(:,:); % Bloque diagonal
      H8x8(6*N_base+1:7*N_base,7*N_base+1:8*N_base) = 0;
      % Octava fila
      H8x8(7*N_base+1:8*N_base,0*N_base+1:1*N_base) = sqrt(1/3)*PP;
      H8x8(7*N_base+1:8*N_base,1*N_base+1:2*N_base) = -sqrt(1/3)*Pz;
      H8x8(7*N_base+1:8*N_base,2*N_base+1:3*N_base) = sqrt(2)*C_d';
      H8x8(7*N_base+1:8*N_base,3*N_base+1:4*N_base) = -sqrt(3/2)*B_d';
      H8x8(7*N_base+1:8*N_base,4*N_base+1:5*N_base) = sqrt(2)*A_d';
      H8x8(7*N_base+1:8*N_base,5*N_base+1:6*N_base) = sqrt(1/2)*B_d;
      H8x8(7*N_base+1:8*N_base,6*N_base+1:7*N_base) = 0;
      H8x8(7*N_base+1:8*N_base,7*N_base+1:8*N_base) = E7PH(:,:); % Bloque diagonal

      ind_filas = (m1+m_angular)*8*N_base+1:(m1+m_angular+1)*8*N_base;
      ind_colum = (m2+m_angular)*8*N_base+1:(m2+m_angular+1)*8*N_base;
      HT(ind_filas,ind_colum) = H8x8(:,:);

    end
  end

  % % Here I force the Hamiltonian to be hermitian
  % HT = triu(HT); % Upper triangular part including diagonal
  % HT = HT + conj(transpose(triu(HT, 1))); % Sum the conjugate transpose of the upper triangular part without the diagonal
  % HT = triu(HT) + triu(HT, 1)';

  for i = 1:8*(2*m_angular+1)
    ST((i-1)*N_base+1:i*N_base,(i-1)*N_base+1:i*N_base) = S(:,:);
  end

  % ST = triu(ST) + triu(ST, 1)';

  auval = eigs(HT, ST, size(HT, 1)-2);
  auval = 1e3*eV*sort(real(auval));

  fprintf(file, '%17.9e  ', k_z/a0);
  for ind_auval = 1:size(auval)
    fprintf(file, '%17.9e  ', auval(ind_auval));
  end
  fprintf(file, '\n');

end
fclose(file);

% quit

% plot(ones(size(auval)), auval,'o')
