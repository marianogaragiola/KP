function v = function_v(h, P0, e, k);

  x(:) = e(:,1)-i*e(:,2);
  x = x';

  v = P0/(sqrt(6)*h)*x*k';

end
