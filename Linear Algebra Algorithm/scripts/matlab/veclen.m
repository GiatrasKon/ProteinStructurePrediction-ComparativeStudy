function vecl = veclen(v,u)
 
N = length(v);
vecl = 0;
 
for i=1:N vecl = vecl + (v(i) - u(i))^2; end
 
vecl = sqrt(vecl);
