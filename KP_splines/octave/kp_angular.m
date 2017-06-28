% Codigo que calcula el hamiltoniano kp segun
% J. Phys.: Condens. Matter 26 (2014) 095501
% pero solo tiene en cuenta el momento angular,
% no mete el campo.
%
% Usa B-splines como funciones base

clear all
close all

format long
pkg load nurbs

a0 = 0.0529177210; eV = 27.21138564; c_light = 137.035999074492;
alpha = 658.4092645439;



% Parametros para los B-splines
num_intervalos   = 8; % numero de subintervalos
N_cuad           = 100;
num_break_points = num_intervalos + 1; % numero de break points contando los extremos
kord             = 5; % orden de los bsplines, el grado de los polinomios es kord-1
regularity       = kord-2; % repeticion de los knots, simplemente degenerados por eso kord-2
xMax             = 20/a0; % extremo superior del intervalo
xMin             = 0; % extremo inferior del intervalo
N_splines        = num_intervalos + kord - 1;
N_base           = N_splines - 1;
num_knots        = N_splines + kord;

% Parametros fisicos.
m0 = 1.0; % electron mass in atomic units
hbar = 1.0; % Planck constant in atomic units.
m_angular = 4; % Momento angular total
k_z = 0*a0; % Componente del vector k en el eje z

N_dim_H = 8*(2*m_angular+1)*N_base;

% Parametros del material (ahora uso de GaAs)
Eg = 1.519/eV; % atomic units
Ep = 25.7/eV; %28.81/eV; % atomic units
delta = 0.33/eV; %0.341/eV; % atomic units
E1 = Eg; % Bottom of the conduction band
E2 = 0.0; % Top of the valence band
P = sqrt(0.5*Ep*hbar^2/m0);
% gamma = [7.05 2.35 3];%
% gamma = [6.85 2.1 2.9];
gamma = [1.8251, -0.2625, 0.3875]; % Son los \tilde{gamma} del paper de Peeters
% gamma = [gamma(1)-Ep/(3*Eg + delta), gamma(2)-0.5*Ep/(3*Eg + delta), gamma(3)-0.5*Ep/(3*Eg + delta)];


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculos las matrices que necesito para los operadores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primero calculo la matriz de solapamiento
S = bs*(repmat(w.*x, 1, N_base).*bs');

% Caculo la matriz de derivadas
T = dbs*(repmat(w.*x, 1, N_base).*dbs');

% Calculo la matriz de rho B_{n1} dB_{n2}/dr
D = bs*(repmat(w.*x, 1, N_base).*dbs');

% Calculo la matriz de B_{n1} dB_{n2}/dr
BdB = bs*(repmat(w, 1, N_base).*dbs');

% Calculo la matriz de rho B_{n1} dB_{n2}/dr
D_rho = bs*(repmat(w.*x.^2,1, N_base).*dbs');

% Matriz del operador rho
rho1 = bs*(repmat(w.*x.^2, 1, N_base).*bs');

% Matriz del operador rho^2
rho2 = bs*(repmat(w.*x.^3, 1, N_base).*bs');

% Matriz del operador 1/rho^2
rho_2 = bs*(repmat(w./x, 1, N_base).*bs');

% Matriz del operador 1/rho
rho_1 = bs*(repmat(w, 1, N_base).*bs');

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
    E7PH = kronDel(m1,m2)*((E2 - delta)*S - 0.5*hbar^2/m0*gamma(1)*(kpkm + k_z^2*S));

    % Calculo el bloque E_{8+}^H
    E8PH = kronDel(m1,m2)*(E2*S - 0.5*hbar^2/m0*((gamma(1) + gamma(2))*kpkm + ...
                                                 (gamma(1) - 2*gamma(2))*k_z^2*S));

    % Calculo el bloque E_{8-}^H
    E8PL = kronDel(m1,m2)*(E2*S - 0.5*hbar^2/m0*((gamma(1) - gamma(2))*kpkm + ...
                                                 (gamma(1) + 2*gamma(2))*k_z^2*S));

    % Calculo el bloque A
    A = kronDel(m1,m2)*(0.5*hbar^2*gamma(2)/m0*(-kpkm - 2.0*k_z^2*S));

    % Calculo el bloque de P^+
    PP = -complex(0,P*kronDel(m2-m1+1,0)*(D - m2*rho_1));

    % Calculo el bloque de P^-
    PM = -complex(0,P*kronDel(m2-m1-1,0)*(D + m2*rho_1));

    % Calculo el bloque de P^z
    Pz = kronDel(m1,m2)*P*k_z*S;

    % Calculo el bloque B
    B = sqrt(3)*hbar^2*gamma(3)/m0*k_z*kronDel(m2-m1-1,0)*(D + m2*rho_1);

    % Calculo el bloque C
    C = 0.25*sqrt(3)*hbar^2/m0*((gamma(2)+gamma(3))*kmkm + (gamma(2)+gamma(3))*kpkp);
    % C = 0.25*sqrt(3)*hbar^2/m0*((gamma(2)+gamma(3))*kmkm + (gamma(2)-gamma(3))*kpkp);

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
    H8x8(2*N_base+1:3*N_base,6*N_base+1:7*N_base) = sqrt(1/2)*B;
    H8x8(2*N_base+1:3*N_base,7*N_base+1:8*N_base) = sqrt(2)*C;
    % Cuarta fila
    H8x8(3*N_base+1:4*N_base,0*N_base+1:1*N_base) = sqrt(2/3)*Pz;
    H8x8(3*N_base+1:4*N_base,1*N_base+1:2*N_base) = -sqrt(1/6)*PM;
    H8x8(3*N_base+1:4*N_base,2*N_base+1:3*N_base) = B';
    H8x8(3*N_base+1:4*N_base,3*N_base+1:4*N_base) = E8PL(:,:); % Bloque diagonal
    H8x8(3*N_base+1:4*N_base,4*N_base+1:5*N_base) = 0;
    H8x8(3*N_base+1:4*N_base,5*N_base+1:6*N_base) = C;
    H8x8(3*N_base+1:4*N_base,6*N_base+1:7*N_base) = -sqrt(2)*A;
    H8x8(3*N_base+1:4*N_base,7*N_base+1:8*N_base) = -sqrt(3/2)*B;
    % Quinta fila
    H8x8(4*N_base+1:5*N_base,0*N_base+1:1*N_base) = sqrt(1/6)*PP;
    H8x8(4*N_base+1:5*N_base,1*N_base+1:2*N_base) = sqrt(2/3)*Pz;
    H8x8(4*N_base+1:5*N_base,2*N_base+1:3*N_base) = C';
    H8x8(4*N_base+1:5*N_base,3*N_base+1:4*N_base) = 0;
    H8x8(4*N_base+1:5*N_base,4*N_base+1:5*N_base) = E8PL(:,:); % Bloque diagonal
    H8x8(4*N_base+1:5*N_base,5*N_base+1:6*N_base) = -B;
    H8x8(4*N_base+1:5*N_base,6*N_base+1:7*N_base) = -sqrt(3/2)*B';
    H8x8(4*N_base+1:5*N_base,7*N_base+1:8*N_base) = sqrt(2)*A;
    % Sexta fila
    H8x8(5*N_base+1:6*N_base,0*N_base+1:1*N_base) = 0;
    H8x8(5*N_base+1:6*N_base,1*N_base+1:2*N_base) = sqrt(1/2)*PP;
    H8x8(5*N_base+1:6*N_base,2*N_base+1:3*N_base) = 0;
    H8x8(5*N_base+1:6*N_base,3*N_base+1:4*N_base) = C';
    H8x8(5*N_base+1:6*N_base,4*N_base+1:5*N_base) = -B';
    H8x8(5*N_base+1:6*N_base,5*N_base+1:6*N_base) = E8PH(:,:); % Bloque diagonal
    H8x8(5*N_base+1:6*N_base,6*N_base+1:7*N_base) = -sqrt(2)*C';
    H8x8(5*N_base+1:6*N_base,7*N_base+1:8*N_base) = sqrt(1/2)*B';
    % Septima fila
    H8x8(6*N_base+1:7*N_base,0*N_base+1:1*N_base) = sqrt(1/3)*Pz;
    H8x8(6*N_base+1:7*N_base,1*N_base+1:2*N_base) = sqrt(1/3)*PM;
    H8x8(6*N_base+1:7*N_base,2*N_base+1:3*N_base) = sqrt(1/2)*B';
    H8x8(6*N_base+1:7*N_base,3*N_base+1:4*N_base) = -sqrt(2)*A';
    H8x8(6*N_base+1:7*N_base,4*N_base+1:5*N_base) = -sqrt(3/2)*B;
    H8x8(6*N_base+1:7*N_base,5*N_base+1:6*N_base) = -sqrt(2)*C;
    H8x8(6*N_base+1:7*N_base,6*N_base+1:7*N_base) = E7PH(:,:); % Bloque diagonal
    H8x8(6*N_base+1:7*N_base,7*N_base+1:8*N_base) = 0;
    % Octava fila
    H8x8(7*N_base+1:8*N_base,0*N_base+1:1*N_base) = sqrt(1/3)*PP;
    H8x8(7*N_base+1:8*N_base,1*N_base+1:2*N_base) = -sqrt(1/3)*Pz;
    H8x8(7*N_base+1:8*N_base,2*N_base+1:3*N_base) = sqrt(2)*C';
    H8x8(7*N_base+1:8*N_base,3*N_base+1:4*N_base) = -sqrt(3/2)*B';
    H8x8(7*N_base+1:8*N_base,4*N_base+1:5*N_base) = sqrt(2)*A';
    H8x8(7*N_base+1:8*N_base,5*N_base+1:6*N_base) = sqrt(1/2)*B;
    H8x8(7*N_base+1:8*N_base,6*N_base+1:7*N_base) = 0;
    H8x8(7*N_base+1:8*N_base,7*N_base+1:8*N_base) = E7PH(:,:); % Bloque diagonal

    ind_filas = (m1+m_angular)*8*N_base+1:(m1+m_angular+1)*8*N_base;
    ind_colum = (m2+m_angular)*8*N_base+1:(m2+m_angular+1)*8*N_base;
    HT(ind_filas,ind_colum) = H8x8(:,:);

  end
end

% Here I force the Hamiltonian to be hermitian
% HT = triu(HT); % Upper triangular part including diagonal
% HT = HT + conj(transpose(triu(HT, 1))); % Sum the conjugate transpose of the upper triangular part without the diagonal

for i = 1:8*(2*m_angular+1)
  ST((i-1)*N_base+1:i*N_base,(i-1)*N_base+1:i*N_base) = S(:,:);
end

% HT = triu(HT) + triu(HT, 1)';
% ST = triu(ST) + triu(ST, 1)';
%
% figure, imagesc(full(real(HT))), colorbar
% figure, imagesc(full(imag(HT))), colorbar

% auval = eigs(HT, ST, size(HT, 1)-1);
auval = eig(HT, ST);
auval = eV*sort(real(auval));

auval = [ones(size(auval)), auval];

save('-ascii', 'auval.dat', 'auval')

% plot(ones(size(auval)), auval,'o')
