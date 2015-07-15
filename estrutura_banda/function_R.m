function R = function_R(h, m0, E_c, E_v, P0, gamma, k);

  R = sqrt(3)*h^2/(2.0*m0)*(gamma(2)*(-k(1)^2+k(2)^2)+i*2.0*gamma(3)*k(1)*k(2));

end
