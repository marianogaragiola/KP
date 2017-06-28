function u = function_v(h, P0, e, k);

  x(:) = e(:,3);
  x = x';

  u = P0/(sqrt(3)*h)*x*k';

end
