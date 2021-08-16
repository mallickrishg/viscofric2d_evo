function Yp = odeViscoelastic_a(~,Y,ss,shz,evl)
% ODE function for faults and shear zones evolution
% Y = [slip,theta,V... - flt
%     [s12,s13,e12,e13] - shz
% We express dY/dt as f(Y) and using RK-4th order to integrate with
% adaptive timesteps
% Rishav Mallick, 2019, EOS

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

if shz.power==1    
    e12dot = (s12)./shz.eta;
    e13dot = (s13)./shz.eta;
elseif shz.power > 1 && shz.combo ~=1
    tau = sqrt(s12.^2 + s13.^2);
    kdis = shz.eta.*(tau.^(shz.power-1));
    e12dot = s12.*kdis;
    e13dot = s13.*kdis;
elseif shz.power > 1 && shz.combo == 1
    % combined flow law - diffusion creep + dislocation creep 
    % (e_ij = A_diff*s_ij + A_dis*(tau^n-1)*s_ij  = A*?_ij*(1 + ?*(|?|^n-1))
    tau = sqrt(s12.^2 + s13.^2);
    kdis = shz.alpha.*shz.eta.*(tau.^(shz.power-1));
    e12dot = s12.*(kdis + shz.eta);
    e13dot = s13.*(kdis + shz.eta);
else
    error('Not a valid rheology - change power/combo')
end

% strain rates
Yp(ss.M*ss.dgf+3:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = e12dot;
Yp(ss.M*ss.dgf+4:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = e13dot;

%% INTERACTIONS
% Acceleration (rate of log(V/Vo))
kv = evl.K*Veff + evl.tau0 + evl.lk1212*e12dot + evl.lk1312*e13dot;
Yp(3:ss.dgf:ss.M*ss.dgf) = (kv - ss.b.*ss.sigma.*dth)./(ss.a.*ss.sigma + damping.*V);

% Stressing Rates (careful with cross-terms - seems like they make the
% numerical solution mess up)
s12dot = 1.*evl.sigma120 + 1.*evl.l1212*e12dot + 1.*evl.l1312*e13dot + evl.kl12*Veff;
s13dot = 1.*evl.sigma130 + 1.*evl.l1213*e12dot + 1.*evl.l1313*e13dot + evl.kl13*Veff;

Yp(ss.M*ss.dgf+1:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = s12dot;
Yp(ss.M*ss.dgf+2:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf) = s13dot;
end


