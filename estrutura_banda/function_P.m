function P = function_P(h, m0, E_c, E_v, P0, gamma, k);

  P = -E_v + gamma(1)*h^2/(2.0*m0)*norm(k)^2;

end
