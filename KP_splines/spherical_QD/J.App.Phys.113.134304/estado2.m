n_estate = 150;

figure,
for i = 500:800

norm = auvec(1:105,i)'*S*auvec(1:105,i);
norm = sqrt(norm);

coef = [auvec(1:105,i); 0]';
phi = bspeval(kord-1, coef, knots, x)/norm;
plot(a0*x, x.^2.*phi'.^2)
disp(i)
pause

end
