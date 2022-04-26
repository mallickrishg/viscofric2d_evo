% Simulates viscoelastic earthquake cycles coupling slip on 
% faults and viscoelastic strain
% beneath the brittle-ductile transition in 2D antiplane (for a single fault).
% Rishav Mallick, 2021, Earth Observatory of Singapore


clear
% close all

addpath ~/Dropbox/scripts/topotoolbox/colormaps/
addpath meshing
addpath odefunction
addpath functions

% Rigidity (MPa)
G = 30e3;

% fault patches
y2i = 0;
M = 40;
y3i = 0;
% base of VW locking (needs to be less than Plate thickness)
basevw = 15e3;


% shear zones
power = 1;% strain rate = A*stress^(power)
burger = 1;% on-1/off-0
% rheological coefficient
etaval = 1e20;%Maxwell viscosity in Pa-s; for power>1, this is A^{-1}

alphaval = 1/100;% for Kelvin element (eta_K = alpha*etaM)

if burger==1
    power = 1;
    etaM = etaval;
    etaK = etaval*alphaval;
end

% shear zone meshing parameters
Nx = 0.2;
Transition = 20e3;% Plate thickness/depth edge of fault domain and beginning of shear zone
% domain size
x3_scale = 30e3;
x2_scale = 300e3;
scf = 1.0; % varies the meshing (leave as is)
%% Create faults and shear zones

% Trimesh shear zone
[ss,shz,~] = create_flt_trishz_src(y2i,y3i,Transition,M,Nx,scf,x3_scale,x2_scale,Transition);
disp(['Number of shear zone elements = ' num2str(shz.N)])
% plate velocity
Vpl = 1e-9; %(m/s)

%% impose earthquake parameters
Teq = 100.*3.15e7;% earthquake every Teq years
ncycles = 10; % number of earthquake cycles

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

% Plot fault and shear zone geometry
figure(2000),clf
subplot(211)
plotpatch(ss,ss.y3c./1e3,1)
hold on,
plotshz(shz,x3c./1e3,1)
axis tight equal, box on
xlabel('x_2 (km)'), ylabel('x_3 (km)');
cb=colorbar;cb.Label.String = 'Depth (km)';cb.Direction='reverse';
set(gca,'YDir','reverse')

subplot(212)
plotpatch(ss,ss.y3c./1e3,1),hold on,
plotshz(shz,x3c./1e3,1)
axis tight equal, box on
xlim([0 x2_scale./1e3])
xlim([0 200])
xlabel('x_2 (km)'), ylabel('x_3 (km)');
cb=colorbar;cb.Label.String = 'Depth (km)';cb.Direction='reverse';
set(gca,'YDir','reverse')
% return
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
% velocity-weakening region is pinned and periodic slip is prescribed
vw = ss.y3c>=0e3 & ss.y3c<=basevw;
ss.b(vw) = ss.a(vw) + 5e-3; % arbitrary number - just make a-b<0
% Fault patches - use appropriate stress conditions on the fault (and taper)
taumax = 3;
eqslip = compute_earthquakeslip(ss,evl,Teq*Vpl,taumax);
taueq = evl.K*eqslip;
taueq(taueq<0) = 0;

% static friction coefficient
ss.mu0 = 0.6*ones(ss.M,1);% this does not matter

% plate velocity (m/s)
ss.Vpl = Vpl.*ones(ss.M,1);

% reference slip rate (m/s)
ss.Vo = 1e-6*ones(ss.M,1);

% shear wave speed (m/s)
ss.Vs = 3e3*ones(ss.M,1);

% Degrees of Freedom [slip, log(v/vo)]
ss.dgf = 2;

% discretization criteria
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
% compute long-term loading rate
[e12pl,e13pl] = edotlongtermsol(x2c./(max(shz.A(:,2))-min(shz.A(:,2))),...
    (x3c-min(shz.A(:,2)))./(max(shz.A(:,2))-min(shz.A(:,2))),...
    power,3000,20);
shz.e12pl = e12pl.*Vpl*0.5/(max(shz.A(:,2))-min(shz.A(:,2)));
shz.e13pl = e13pl.*Vpl*0.5/(max(shz.A(:,2))-min(shz.A(:,2)));

% Effective Viscosity (MPa s)
shz.power = power;

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

figure(1),clf
toplot = zeros(shz.N,1);
toplot(shz.e12pl<0) = 1;
plotshz(shz,toplot,1)
axis tight equal
% caxis([1e-17 1e-13])
% caxis([-1 1].*1e-15)
% xlim([0 1].*20)
set(gca,'ColorScale','log','YDir','reverse')
% return
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%         N U M E R I C A L   S O L U T I O N           %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Initialize State Vector
Y0=zeros(ss.M*ss.dgf+shz.N*shz.dgf,1);

Y0(1:ss.dgf:ss.M*ss.dgf) = zeros(ss.M,1);% slip
Y0(2:ss.dgf:ss.M*ss.dgf) = log(ss.Vpl*0.99./ss.Vo); % v

% Shear zones
Y0(ss.M*ss.dgf+1:shz.dgf:end) = 1e-9 + zeros(shz.N,1); %stress12
Y0(ss.M*ss.dgf+2:shz.dgf:end) = zeros(shz.N,1); %stress13
Y0(ss.M*ss.dgf+3:shz.dgf:end) = zeros(shz.N,1); %strain12t
Y0(ss.M*ss.dgf+4:shz.dgf:end) = zeros(shz.N,1); %strain13t
if burger==1
    Y0(ss.M*ss.dgf+5:shz.dgf:end) = zeros(shz.N,1); %strain12k
    Y0(ss.M*ss.dgf+6:shz.dgf:end) = zeros(shz.N,1); %strain13k
end

%% Simulation 
disp('Solving ODE')

tic
% initialize the function handle with
% set constitutive parameters
if burger==0
    yp=@(t,y) odeimposedViscoelastic_a(t,y,ss,shz,evl);
else
    yp=@(t,y) odeimposedViscoelastic_burger(t,y,ss,shz,evl);
end
tic
% Solve the system
options=odeset('Refine',1,'AbsTol',1e-9,'RelTol',1e-9,'InitialStep',1e-6,'MaxStep',3e8); 
for i = 1:ncycles
    if i==1
        [t,Y]=ode45(yp,[0 Teq],Y0,options);
    else
        % provide new initial conditions (loaded by earthquake)
        Y0 = Y(end,:);
        Y0(1:ss.dgf:ss.M*ss.dgf) = Y0(1:ss.dgf:ss.M*ss.dgf)' + eqslip;        
        Y0(2:ss.dgf:ss.M*ss.dgf) = Y0(2:ss.dgf:ss.M*ss.dgf)' + taueq./(ss.a-ss.b)./ss.sigma;
        
        Y0(ss.M*ss.dgf+1:shz.dgf:end) = Y0(ss.M*ss.dgf+1:shz.dgf:end)' + evl.kl12*eqslip; %stress12
        Y0(ss.M*ss.dgf+2:shz.dgf:end) = Y0(ss.M*ss.dgf+2:shz.dgf:end)' + evl.kl13*eqslip; %stress13
        
        [tmod,Ymod]=ode45(yp,[0 Teq],Y0,options);
        t = tmod;
        Y = Ymod;
    end
    disp(['Cycle ' num2str(i) ' in progress'])
end
toc

%% results
V=repmat(ss.Vo',size(Y,1),1).*exp(Y(:,2:ss.dgf:ss.M*ss.dgf));
V(:,vw) = 0;
slip = Y(:,1:ss.dgf:ss.M*ss.dgf);

s12 = Y(:,ss.M*ss.dgf+1:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
s13 = Y(:,ss.M*ss.dgf+2:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);

if burger==0
    e12 = Y(:,ss.M*ss.dgf+3:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
    e13 = Y(:,ss.M*ss.dgf+4:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
    if power==1
        e12d = s12./repmat(shz.eta',size(Y,1),1);
        e13d = s13./repmat(shz.eta',size(Y,1),1);
        viscosity_eff = repmat(shz.eta',length(t),1).*1e6;
    elseif power>1
        j2 = sqrt(s12.^2 + s13.^2);
        e12d = repmat(shz.eta',size(Y,1),1).*s12.*(j2.^(power-1));
        e13d = repmat(shz.eta',size(Y,1),1).*s13.*(j2.^(power-1));
        viscosity_eff = 1e6./(repmat(shz.eta',length(t),1).*j2.^(shz.power-1));
    else
        disp('Invalid Rheology')
    end
else
    e12K = Y(:,ss.M*ss.dgf+5:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
    e13K = Y(:,ss.M*ss.dgf+6:shz.dgf:ss.M*ss.dgf+shz.N*shz.dgf);
    
    %e12d_K = (s12 + (evl.l1212*e12K' + evl.l1312*e13K')')./repmat((shz.etaK)',size(Y,1),1);
    %e13d_K = (s13 + (evl.l1213*e12K' + evl.l1313*e13K')')./repmat((shz.etaK)',size(Y,1),1);
    e12d_K = (s12 - G.*e12K)./repmat((shz.etaK)',size(Y,1),1);
    e13d_K = (s13 - G.*e13K)./repmat((shz.etaK)',size(Y,1),1);
    
    e12d = s12./repmat(shz.etaM',size(Y,1),1) + e12d_K;
    e13d = s13./repmat(shz.etaM',size(Y,1),1) + e13d_K;
end

toc
disp('ODE Solution')
%% PLOT snapshot transects at various t/T_eq

% surface velocties
if isfield(shz,'L')
    deepx3 = max(shz.xx3c(:));
elseif isfield(shz,'tri')
    deepx3 = max(max([shz.A(:,2) shz.B(:,2) shz.C(:,2)]));
end
% deepx3 = 20e3;

% surface dispalcements/velocities
ox = linspace(0.01e3,200e3,100)';
obs = [ox,zeros(length(ox),1)];
Gd = compute_displacementkernels(obs,ss,shz);

vsurf_f = (Gd.kd*V')';
vsurf_deep = (repmat(Vpl/pi.*atan2(ox,deepx3),1,length(t)))';
vsurf_v = (Gd.l12d*e12d' +  Gd.l13d*e13d')'; % no deep loading
% vsurf_v = (Gd.l12d*(e12d'-shz.e12pl) +  Gd.l13d*(e13d'-shz.e13pl))'; % no deep loading
vsurf = vsurf_v + vsurf_deep;


figure(10),clf
set(gcf,'Position',[0 0.5 2 2].*500)
% tplotvec = Tevent + [1/50,1/20,.1,.2,.3,.4,.5,.6,.7,.8,.9].*Teq;
tplotvec = [[0,1,7,30,60,90,180,365,2*365,5*365].*86400,[0.1:0.1:1].*Teq];

cspec = [parula(10);cool(length(tplotvec)-6)];

for i = 1:length(tplotvec)
    index = find(abs(t-tplotvec(i))==min(abs(t-tplotvec(i))),1);
    plot(ox./20e3,vsurf(index,:)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
end
plot(ox./20e3,1/pi.*atan2(ox,20e3),'k-','LineWidth',2)
% plot(ox./1e3,1/pi.*atan2(ox,30e3),'k--','LineWidth',1)
% plot(ox./20e3,1/pi.*atan2(ox,40e3),'k--','LineWidth',1)
% plot(ox./20e3,1/pi.*atan2(ox,60e3),'k--','LineWidth',1)
xlabel('x_2/D'), ylabel('V/V_{pl}')
axis tight, grid on
ylim([1e-2 10])
set(gca,'FontSize',20,'LineWidth',2,'YScale','log')
if burger==0
    title(['Power n = ' num2str(power)])
else
    title('Burgers')
end

%% solve for slip rate and locking depth for final vsurf

% locking depth function
func = @(beta,obsx) beta(1)./pi.*atan2(obsx,beta(2));
func2 = @(beta,obsx,data) data - beta(1)./pi.*atan2(obsx,beta(2));

beta0 = [1,1];
[beta1,res,~,covb,~] = nlinfit(ox./20e3,vsurf(end,:)'./Vpl,func,beta0);
beta2 = lsqnonlin(@(beta) func2(beta,ox./20e3,vsurf(end,:)'./Vpl),beta0,[0.9,0],[1.1,100]);
% betaCI = nlparci(beta1,res,'covar',covb);

figure(5),clf
plot(ox./20e3,vsurf(end,:)'./Vpl,'p','LineWidth',1), hold on
plot(ox./20e3,func(beta1,ox./20e3),'-','LineWidth',2,'Color',rgb('gray'))
plot(ox./20e3,func(beta2,ox./20e3),'r-','LineWidth',2)
plot(ox./20e3,func(beta0,ox./20e3),'k-','LineWidth',2)
% plot(ox./20e3,func([0.9,2],ox./20e3),'k-.','LineWidth',2)
axis tight, grid off
xlabel('x_2/D'), ylabel('v/V_{pl}')
ylim([0 0.5])
legend('interseismic velocity',['V_{pl} = ' num2str(beta1(1),'%.1f') ', D = ' num2str(beta1(2),'%.1f') ' (unbounded)']...
    ,['V_{pl} = ' num2str(beta2(1),'%.1f') ', D = ' num2str(beta2(2),'%.1f') ' (bounded)'],...
    'Elastic Plate',...
    'location','best')
set(gca,'FontSize',15,'LineWidth',1.5)

% T = table(ox./20e3,vsurf(end,:)'./Vpl);
% writetable(T,'Figures/burgers_interseismic.dat');
%% plot snapshots

tplotvec = [0.02,0.95].*Teq;

figure(11),clf
set(gcf,'Position',[0 0.5 3 1.7].*500)

for i = 1:length(tplotvec)
    index = find(abs(t-tplotvec(i))==min(abs(t-tplotvec(i))),1);
    subplot(2,1,i)


%     toplot = sqrt((e12d(index,:)-shz.e12pl').^2 + (e13d(index,:)-shz.e13pl').^2)';
%     toplot = abs(sqrt(e12d(index,:)'.^2 + e13d(index,:)'.^2) - sqrt(shz.e12pl.^2 + shz.e13pl.^2));
%     toplot = abs(e12d(index,:)' - shz.e12pl);
%     toplot = sqrt(shz.e12pl.^2 + shz.e13pl.^2);
    toplot = sqrt(e12d(index,:).^2 + e13d(index,:).^2)';
    
    F = scatteredInterpolant(x2c,x3c,(toplot),'natural','none');
    x2g = linspace(-100e3,100e3,400)';
    %x3g = linspace(min(shz.A(:,2)),max(shz.A(:,2)),200)';
    x3g = linspace(Transition,Transition+x3_scale,200)';
    [X2g,X3g] = meshgrid(x2g,x3g);
    toplotint = F(X2g(:),X3g(:));
    
%     imagesc(x2g./1e3,x3g./1e3,reshape(toplotint,length(x3g),length(x2g))), shading flat, hold on
    plotshz(shz,toplot,1), shading flat, box on, hold on
    contour(x2g./1e3,x3g./1e3,reshape(toplotint,length(x3g),length(x2g)),...
        logspace(log10(Vpl*1e-6),log10(Vpl*1e-3),20),'w-','LineWidth',.1)
    set(gca,'YDir','reverse','FontSize',20,'Color','none','LineWidth',2,'ColorScale','log')
    xlabel('x_2 (km)'),ylabel('x_3 (km)'),
    caxis(Vpl.*[1e-6 1e-4])
    axis tight equal
    xlim([-1 1].*100),    
    ylim([Transition Transition+x3_scale]./1e3)
    colormap(ttscm('oslo',500))
    cb=colorbar;cb.Label.String = 'strain rate';
    title(['$\frac{\Delta t}{T_{eq}} = $' num2str(round((t(index)))./Teq)],'interpreter','latex','Fontsize',25)
end

