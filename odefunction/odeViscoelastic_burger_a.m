function Yp = odeViscoelastic_burger_a(~,Y,ss,shz,evl)
% ODE function for faults and shear zones evolution (Burger's rheology)
% Y = [slip,theta,V... - flt
%     [s12,s13,e12t,e13t,e12k,e13k] - shz
% We express dY/dt as f(Y) and using RK-4th order to integrate with
% adaptive timesteps
% Rishav Mallick, 2021, EOS

%% FAULTS
G = 30e3; % shear modulus
damping = G./ss.Vs/2;

% State variables
th   = Y(2:ss.dgf:ss.M*ss.dgf);

% Slip velocities 
V    = ss.Vo.*exp(Y(3:ss.dgf:ss.M*ss.dgf));
if isfield(ss,'Lf')
    Veff = V./2./ss.Lf;
else
    Veff = V;
end

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1:ss.dgf:ss.M*ss.dgf) = V;

% Rate of state (rate of log(theta/theta0))
dth=(ss.Vo.*exp(-th) - V)./ss.l;
Yp(2:ss.dgf:ss.M*ss.dgf) = dth;

%% SHEAR ZONES
s12 = Y(ss.M*ss.dgf+1:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
s13 = Y(ss.M*ss.dgf+2:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
e12_K = Y(ss.M*ss.dgf+5:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
e13_K = Y(ss.M*ss.dgf+6:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);

% kill cross terms for kelvin strain rates
e12dot_K = (s12 + (evl.l1212*e12_K + evl.l1312*e13_K))./(shz.alpha.*shz.eta);
e13dot_K = (s13 + (evl.l1213*e12_K + evl.l1313*e13_K))./(shz.alpha.*shz.eta);

e12dot_T = (s12)./shz.eta + e12dot_K;
e13dot_T = (s13)./shz.eta + e13dot_K;

% strain rates
% total
Yp(ss.M*ss.dgf+3:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = e12dot_T;
Yp(ss.M*ss.dgf+4:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = e13dot_T;

%kelvin
Yp(ss.M*ss.dgf+5:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = e12dot_K;
Yp(ss.M*ss.dgf+6:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = e13dot_K;
%% INTERACTIONS
% Acceleration (rate of log(V/Vo))
kv = evl.K*Veff + evl.tau0 + evl.lk1212*e12dot_T + evl.lk1312*e13dot_T;
Yp(3:ss.dgf:ss.M*ss.dgf) = (kv - ss.b.*ss.sigma.*dth)./(ss.a.*ss.sigma + damping.*V);

% Stressing Rates (careful with cross-terms - seems like they make the
% numerical solution mess up)
s12dot = evl.sigma120 + evl.l1212*e12dot_T + evl.l1312*e13dot_T + evl.kl12*Veff;
s13dot = evl.sigma130 + evl.l1213*e12dot_T + evl.l1313*e13dot_T + evl.kl13*Veff;

Yp(ss.M*ss.dgf+1:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = s12dot;
Yp(ss.M*ss.dgf+2:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = s13dot;
end