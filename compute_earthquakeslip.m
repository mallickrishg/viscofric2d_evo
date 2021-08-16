function slip_fullfault = compute_earthquakeslip(ss,evl,maxslip,taumax)

vw = ss.a<ss.b;
vs = ~vw;

slip_fullfault = zeros(ss.M,1);
slip_fullfault(vw) = maxslip;
slipco = slip_fullfault(vw);

Nco = length(slipco);
Naf = ss.M-Nco;

% set bounds on slip amplitude
lb = zeros(Naf,1);
ub = maxslip.*ones(Naf,1);

% stress kernels (sub-sets)
Kaf = evl.K(vs,vs);
Kco = evl.K(vs,vw);

tauco = Kco*slipco;

Ain = [Kaf;-Kaf];
bin = [(taumax.*ones(Naf,1) - tauco);tauco];

slip_af = lsqlin(eye(Naf),zeros(Naf,1),Ain,bin,[],[],lb,ub);
slip_fullfault(vs) = slip_af;

end