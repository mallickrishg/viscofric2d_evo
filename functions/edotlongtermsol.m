function [e12d,e13d] = edotlongtermsol(x,z,n,iter,w)
S12 = zeros(size(x));
S13 = S12;

for m = 1:iter
    S12 = S12 + cosh(m*pi*(1-z)./w./sqrt(n))./cosh(m*pi/w/sqrt(n)).*cos(m*pi*x/w);
    
    S13 = S13 + sinh(m*pi*(1-z)./w./sqrt(n))./cosh(m*pi/w/sqrt(n)).*sin(m*pi*x/w);
end
e12d = 1/(2*w) + 1/w*S12;
e13d = -1/w/sqrt(n)*S13;

end