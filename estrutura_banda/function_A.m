function A = function_A(h, m0, E_c, E_v, P0, gamma, k);

  A = E_c + h^2/(2.0*m0)*norm(k)^2;

end
