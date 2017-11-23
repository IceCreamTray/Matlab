%%
% Boyfriend contribution

lotsa = -10 : 0.01 : 10;
lo = 16 * (power(sin(lotsa), 3));
ve = (13 * cos(lotsa)) - (5 * cos(2 * lotsa)) - (2 * cos(3 * lotsa)) - (cos(4 * lotsa));
	
figure
plot(lo,ve);
title('929');
