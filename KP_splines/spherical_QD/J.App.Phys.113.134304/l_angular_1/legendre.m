function [y] = legendre(x, n)
  %legendre recursively computes the value of x for the Legendre  polynomial n
  
  if n == 0 
      %y = ones(1,length(x));
      y = ones(size(x));
  elseif n == 1
      y = x;
  else
      y = ((2*(n-1)+1) * x .* legendre(x,n-1) - ((n-1) * legendre(x,n-2)))/(n);
  end

end
