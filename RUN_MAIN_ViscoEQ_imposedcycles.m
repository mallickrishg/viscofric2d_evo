function RUN_MAIN_ViscoEQ_imposedcycles(Nx,M,Transition,burger,power,etaval,Trecur,Vpl)
% Nx - shear zone mesh number
% M - number of fault patches
% Transition - elastic layer thickness (m) generally 20e3
% power, etaval - power and rheological coefficient (in Pa-s)
% Trecur - recurrence time of earthquake in years
% Vpl - plate rate in (m/s) - generally 1e-9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% burger - (0/1) to turn on burger's biviscous rheology
% in the case of burger's rheology, 
% etaval - Maxwell viscosity
% etaval*power - Kelvin viscosity

addpath functions/
addpath odefunction/

% Rigidity (MPa)
G = 30e3;

% fault patches
y2i = 0;
% M = 100;
y3i = 0;

% shear zones
combo = 0;% use diffusion + dislocation creep

% etaval = 1e18;%A_0 = 1/?
% use rheology of the form - e_ij = A_0*sigma_ij + A_0*alpha*(|sigma|^n-1)*sigma_ij = A_0*sigma_ij*(1 + alpha*(|sigma|^n-1))
alphaval = power;% relative weight of diffusion creep and dislocation creep

if burger==1
    power = 1;
    etaM = etaval;
    etaK = etaval*alphaval;
end

% Nx = 50;
% Transition = 20e3;
Transitionshz = Transition+0e3;

% base of VW locking
basevw = 15e3;

x3_scale = 30e3;
x2_scale = 390e3;% 590 km for Tr = 200, 390 for Tr - 50, for power-law - 190 is good enough

scf = 1.0;
%% Create faults and shear zones
% [ss,shz,~] = create_flt_shz_src(y2i,y3i,Transition,M,Nx,Nz,x3_scale,x2_scale,Transitionshz);

% Trimesh shear zone
[ss,shz,~] = create_flt_trishz_src(y2i,y3i,Transition,M,Nx,scf,x3_scale,x2_scale,Transitionshz);

y2i = unique(ss.y2f);

% plate velocity
% Vpl = 1e-9; %(m/s)

%% impose earthquake parameters
Teq = Trecur.*3.15e7;% earthquake every Teq years
ncycles = 10000/Trecur; % number of earthquake cycles

%% Stress Kernels and EVL object

evl = compute_stresskernels(ss,shz);
disp('Loaded Stress Kernels')


if isfield(shz,'L')
    x2c = shz.xx2c(:);
    x3c = shz.xx3c(:);
elseif isfield(shz,'tri')
    x2c = mean([shz.A(:,1),shz.B(:,1),shz.C(:,1)],2);
    x3c = mean([shz.A(:,2),shz.B(:,2),shz.C(:,2)],2);
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% effective confining pressure on fault (MPa)
ss.sigma = 50*ones(ss.M,1);

% frictional parameters
ss.a = 1e-2*ones(ss.M,1);
ss.b = ss.a - 5e-3*ones(ss.M,1);
% velocity-weakening
vw = ss.y3c>=2e3 & ss.y3c<=basevw;
ss.b(vw) = ss.a(vw) + 5e-3;

% plate velocity (m/s)
ss.Vpl = Vpl.*ones(ss.M,1);

% reference slip rate (m/s)
ss.Vo = 1e-6*ones(ss.M,1);

% shear wave speed (m/s)
ss.Vs = 3e3*ones(ss.M,1);

% Degrees of Freedom [slip, log(v/vo)]
ss.dgf = 2;

dx_fault = (max(ss.Wf));
if isfield(shz,'L')
    dx_shearzone = (max(shz.W(:)));
else
    dx_shearzone = min(sqrt((shz.A(:,1)-shz.B(:,1)).^2 + (shz.A(:,2)-shz.B(:,2)).^2));
end

disp(['dx(fault) = ' num2str(dx_fault) 'm'])
disp(['dx(shear) = ' num2str(dx_shearzone) 'm'])

disp('Assigned Frictional Properties')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                   R H E O L O G Y                    %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
[e12pl,e13pl] = edotlongtermsol(x2c./(max(shz.A(:,2))-min(shz.A(:,2))),...
    (x3c-min(shz.A(:,2)))./(max(shz.A(:,2))-min(shz.A(:,2))),...
    power,1000,10);
shz.e12pl = e12pl.*Vpl*0.5/(max(shz.A(:,2))-min(shz.A(:,2)));
shz.e13pl = e13pl.*Vpl*0.5/(max(shz.A(:,2))-min(shz.A(:,2)));

% Effective Viscosity (MPa s)
shz.power = power;
shz.combo = combo;
shz.alpha = alphaval;

if burger == 0
    if power==1
        shz.eta = etaval*1e-6.*ones(shz.N,1);% convert to MPa s
    elseif power>1
        %here eta is actually 1/eta_eff
        shz.eta = 1./(etaval*1e-6).*ones(shz.N,1);
    else
        disp('Unrecognized rheology')
    end
    % degrees of freedom
    shz.dgf = 4;
else % burgers rheology (linear)
    shz.etaM = etaM.*1e-6.*ones(shz.N,1);
    shz.etaK = etaK.*1e-6.*ones(shz.N,1);
    % degrees of freedom
    shz.dgf = 6;
end
disp('Assigned Rheological Properties')


% Loading Stresses
% backslip loading
evl.tau0     = -evl.K*ss.Vpl    - evl.lk1212*shz.e12pl - evl.lk1312*shz.e13pl;
evl.sigma120 = -evl.kl12*ss.Vpl - evl.l1212*shz.e12pl  - evl.l1312*shz.e13pl;
evl.sigma130 = -evl.kl13*ss.Vpl - evl.l1213*shz.e12pl  - evl.l1313*shz.e13pl;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%         N U M E R I C A L   S O L U T I O N           %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% output directory
if burger==0
    odir = ['imposedviscocycles_' num2str(Vpl,'%1.e') '/Nx_' num2str(Nx) '_M_' num2str(M) ...
    '/power_' num2str(power) '_' num2str(etaval) '/Trec_' num2str(Teq./3.15e7)];
else
    odir = ['imposedviscocycles_' num2str(Vpl,'%1.e') '/Nx_' num2str(Nx) '_M_' num2str(M) ...
    '/burger_' num2str(etaM) '_' num2str(etaK) '/Trec_' num2str(Teq./3.15e7)];
end

mkdir(odir)

save([odir '/modelgeometry.mat'],'ss','shz','evl','Teq')

disp(odir)

%% solve ODE
disp('Solving ODE')

% Initialize State Vector
Y0=zeros(ss.M*ss.dgf+shz.N*shz.dgf,1);

% Fault patches
taumax = 3;% maximum shear stress increase in VS region

eqslip = compute_earthquakeslip(ss,evl,Teq*Vpl,taumax);
taueq = evl.K*eqslip;
taueq(taueq<0) = 0;

Y0(1:ss.dgf:ss.M*ss.dgf) = zeros(ss.M,1);% slip
Y0(2:ss.dgf:ss.M*ss.dgf) = log(ss.Vpl*0.99./ss.Vo); % v

% Shear zones
Y0(ss.M*ss.dgf+1:shz.dgf:end) = 1e-9 + zeros(shz.N,1); %stress12
Y0(ss.M*ss.dgf+2:shz.dgf:end) = zeros(shz.N,1); %stress13
Y0(ss.M*ss.dgf+3:shz.dgf:end) = zeros(shz.N,1); %strain12
Y0(ss.M*ss.dgf+4:shz.dgf:end) = zeros(shz.N,1); %strain13

%% Simulation 

tic
% initialize the function handle with
if burger==0
    yp=@(t,y) odeimposedViscoelastic_a(t,y,ss,shz,evl);
else
    yp=@(t,y) odeimposedViscoelastic_burger(t,y,ss,shz,evl);
end
tic
% Solve the system
options=odeset('Refine',1,'AbsTol',1e-6,'RelTol',1e-6,'InitialStep',1e-6,'MaxStep',3e8); 
for i = 1:ncycles
    if i==1
        [t,Y]=ode45(yp,[0 Teq],Y0,options);
        %t=t';Y=Y';
    else
        % provide new initial conditions (loaded by earthquake)
        Y0 = Y(end,:);
        Y0(1:ss.dgf:ss.M*ss.dgf) = Y0(1:ss.dgf:ss.M*ss.dgf)' + eqslip;        
        Y0(2:ss.dgf:ss.M*ss.dgf) = Y0(2:ss.dgf:ss.M*ss.dgf)' + taueq./(ss.a-ss.b)./ss.sigma;
        
        Y0(ss.M*ss.dgf+1:shz.dgf:end) = Y0(ss.M*ss.dgf+1:shz.dgf:end)' + evl.kl12*eqslip; %stress12
        Y0(ss.M*ss.dgf+2:shz.dgf:end) = Y0(ss.M*ss.dgf+2:shz.dgf:end)' + evl.kl13*eqslip; %stress13
        
        [t,Y]=ode45(yp,[0 Teq],Y0,options);
        %t=[t;(i-1)*Teq+tmod'];Y=[Y;Ymod'];
        
    end
    disp(['Cycle ' num2str(i) ' in progress'])
end
toc

%% results
V=repmat(ss.Vo',size(Y,1),1).*exp(Y(:,2:ss.dgf:ss.M*ss.dgf));
V(:,vw) = nan;
slip = Y(:,1:ss.dgf:ss.M*ss.dgf);

s12 = Y(:,ss.M*ss.dgf+1:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
s13 = Y(:,ss.M*ss.dgf+2:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);

if burger==0
    e12 = Y(:,ss.M*ss.dgf+3:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
    e13 = Y(:,ss.M*ss.dgf+4:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
    if power==1
        e12d = s12./repmat(shz.eta',size(Y,1),1);
        e13d = s13./repmat(shz.eta',size(Y,1),1);
        %viscosity_eff = repmat(shz.eta',length(t),1).*1e6;
    elseif power>1
        j2 = sqrt(s12.^2 + s13.^2);
        e12d = repmat(shz.eta',size(Y,1),1).*s12.*(j2.^(power-1));
        e13d = repmat(shz.eta',size(Y,1),1).*s13.*(j2.^(power-1));
        %viscosity_eff = 1e6./(repmat(shz.eta',length(t),1).*j2.^(shz.power-1));
    else
        disp('Invalid Rheology')
    end
    save([odir '/modeloutputs_s.mat'],'e12d','e13d','s12','s13','-v7.3')
else
    e12K = Y(:,ss.M*ss.dgf+5:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
    e13K = Y(:,ss.M*ss.dgf+6:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
    
    %e12d_K = (s12 + (evl.l1212*e12K' + evl.l1312*e13K')')./repmat((shz.etaK)',size(Y,1),1);
    %e13d_K = (s13 + (evl.l1213*e12K' + evl.l1313*e13K')')./repmat((shz.etaK)',size(Y,1),1);
    e12d_K = (s12 - G.*e12K)./repmat((shz.etaK)',size(Y,1),1);
    e13d_K = (s13 - G.*e13K)./repmat((shz.etaK)',size(Y,1),1);
    
    e12d = s12./repmat(shz.etaM',size(Y,1),1) + e12d_K;
    e13d = s13./repmat(shz.etaM',size(Y,1),1) + e13d_K;
    save([odir '/modeloutputs_s.mat'],'e12d','e13d','e12d_K','e13d_K','s12','s13','-v7.3')
end

%%%%%%%%%%%%%%%%%%%%%% SAVE OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%
save([odir '/modeloutputs_f.mat'],'t','V','slip','-v7.3')


toc
disp('ODE Solution')

end