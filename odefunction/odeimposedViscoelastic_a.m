function Yp = odeimposedViscoelastic_a(~,Y,ss,shz,evl)
% ODE function for faults and shear zones evolution
% Y = [slip,V... - flt
%     [s12,s13,e12,e13] - shz
% We express dY/dt as f(Y) and using RK-4th order to integrate with
% adaptive timesteps
% Rishav Mallick, 2019, EOS

%% FAULTS
G = 30e3; % shear modulus
damping = G./ss.Vs/2;
VS = ss.a>ss.b;

% Slip velocities 
dummy = ss.Vo.*exp(Y(2:ss.dgf:ss.M*ss.dgf));
V = zeros(ss.M,1);
V(VS) = dummy(VS);

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1:ss.dgf:ss.M*ss.dgf) = V;

%% SHEAR ZONES
s12 = Y(ss.M*ss.dgf+1:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
s13 = Y(ss.M*ss.dgf+2:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);

if shz.power==1    
    e12dot = (s12)./shz.eta;
    e13dot = (s13)./shz.eta;
elseif shz.power > 1
    tau = sqrt(s12.^2 + s13.^2);
    kdis = shz.eta.*(tau.^(shz.power-1));
    e12dot = s12.*kdis;
    e13dot = s13.*kdis;
else
    error('Not a valid rheology - change power/combo')
end

% strain rates
Yp(ss.M*ss.dgf+3:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = e12dot;
Yp(ss.M*ss.dgf+4:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = e13dot;

%% INTERACTIONS
% Acceleration (rate of log(V/Vo))
kv = evl.K*V + evl.tau0 + evl.lk1212*e12dot + evl.lk1312*e13dot;
kv(~VS) = 0;
Yp(2:ss.dgf:ss.M*ss.dgf) = (kv)./((ss.a-ss.b).*ss.sigma + damping.*V);

% Stressing Rates (careful with cross-terms - seems like they make the
% numerical solution mess up)
s12dot = evl.sigma120 + evl.l1212*e12dot + evl.l1312*e13dot + evl.kl12*V;
s13dot = evl.sigma130 + evl.l1213*e12dot + evl.l1313*e13dot + evl.kl13*V;

Yp(ss.M*ss.dgf+1:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = s12dot;
Yp(ss.M*ss.dgf+2:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = s13dot;
end


