n_estate = 650;

figure, hold on

norm = auvec(1:105,n_estate)'*S*auvec(1:105,n_estate);
norm = sqrt(norm);

coef = [auvec(1:105,n_estate); 0]';

phi = bspeval(kord-1, coef, knots, x)/norm;
plot(a0*x, x.^2.*phi'.^2)
