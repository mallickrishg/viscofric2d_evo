function u1d = udotlongtermsol(x,z,n,iter,w)
S1 = zeros(size(x));

for m = 1:iter
    S1 = S1 + cosh(m*pi*(1-z)./w./sqrt(n))./cosh(m*pi/w/sqrt(n)).*sin(m*pi*x/w)./m;    
end
u1d = x./w + 2/pi.*S1;

end