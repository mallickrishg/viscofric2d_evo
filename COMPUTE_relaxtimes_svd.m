% load outputs from simulations and compute relaxation time for the system
% using SVD
% Rishav Mallick, Caltech Seismolab, 2022

clear
addpath functions/
addpath ~/Dropbox/scripts/utils/
addpath ~/Dropbox/scripts/topotoolbox/colormaps/

% specify parameters
Nx = 0.2;
M = 100;
Vpl = 1e-9;
Trecur = 100;
burger = 0;
etaval = 1e18;
powerval = 1;

if burger==0
    odir = ['imposedviscocycles_' num2str(Vpl,'%1.e') '/Nx_' num2str(Nx) '_M_' num2str(M) ...
    '/power_' num2str(powerval) '_' num2str(etaval) '/Trec_' num2str(Trecur)];
else
    odir = ['imposedviscocycles_' num2str(Vpl,'%1.e') '/Nx_' num2str(Nx) '_M_' num2str(M) ...
    '/burger_' num2str(etaval) '_' num2str(etaval*powerval) '/Trec_' num2str(Trecur)];
end

disp(odir)

% load results
load([odir '/modelgeometry.mat'],'ss','shz','evl','Teq')
load([odir '/modeloutputs_f.mat'])
load([odir '/modeloutputs_s.mat'])

%% stress vs displacement and velocity
tindex = t<=3.15e7*5;
V(isnan(V)) = 0;
sigmadot_f = evl.K*(V'-Vpl) + evl.lk1212*(e12d'-shz.e12pl) + evl.lk1312*(e13d'-shz.e13pl);
sigma_f = cumtrapz(t,sigmadot_f');
slip = cumtrapz(t,V);

figure(1),clf
imagesc(sigma_f')
colorbar;
caxis([-1 1]*5)
colormap(ttscm('oleron',50))

figure(2),clf
fid = [10,80,90];
plot(slip(tindex,fid),sigma_f(tindex,fid),'-','Linewidth',2)
axis tight
xlim([0 1])
ylim([-3 0])

figure(3),clf
semilogx(V(tindex,fid),sigma_f(tindex,fid),'o-','Linewidth',2)
axis tight
xlim([0.1 1e5].*Vpl)
ylim([-3 0])

% shear zone
x2c = mean([shz.A(:,1),shz.B(:,1),shz.C(:,1)],2);
x3c = mean([shz.A(:,2),shz.B(:,2),shz.C(:,2)],2);
x20 = 0;x30 = 20e3;
% distance from fault corner (m)
r = sqrt((x2c-x20).^2 + (x3c-x30).^2);
[~,rind] = sort(r);
triarea = 0.5*abs(shz.A(:,1).*(shz.B(:,2) - shz.C(:,2)) + ...
    shz.B(:,1).*(shz.C(:,2)-shz.A(:,2)) + shz.C(:,1).*(shz.A(:,2)-shz.B(:,2)));

figure(4),clf
sid = r>=1e3 &r<=2e3;
xplot = e12d(tindex,sid)*triarea(sid)/sum(triarea(sid));
yplot = s12(tindex,sid)*triarea(sid)/sum(triarea(sid));
plot(xplot,yplot,'.','Linewidth',1)
axis tight, grid on

e12 = cumtrapz(t,e12d);
figure(5),clf
xplot = e12(tindex,sid)*triarea(sid)/sum(triarea(sid));
plot(xplot,yplot,'-','Linewidth',2)
axis tight
xlim([0 1e-3])
ylim([0 10])

return
%% decompose timeseries using SVD
% displacement greens functions
ox = linspace(-200e3,200e3,200)';
obs = [ox,zeros(length(ox),1)];
Gd = compute_displacementkernels(obs,ss,shz);

% surface velocties
vsurf_v = (Gd.l12d*(e12d-shz.e12pl')' +  Gd.l13d*(e13d-shz.e13pl')')'; % no deep loading
% surface velocities
vsurf = vsurf_v;
% surface displacement
usurf= cumtrapz(t,vsurf);

% use SVD on displacements
twindow = 4;% in years
plotindex = t<=twindow*3.15e7;
[U,S,V] = svd(usurf(plotindex,:),'econ');

figure(10),clf
subplot(121)
semilogy(1:length(diag(S)), diag(S)./sum(diag(S)),'x-','Linewidth',2)
axis tight, box on, grid on

subplot(122)
plot(1:length(diag(S)),cumsum(diag(S))./sum(diag(S)),'o-','Linewidth',2)
axis tight, box on, grid on
ylim([0 1])
xlim([0 10])

figure(11),clf
svdnum = 1;
subplot(121)
plot(ox./1e3,S(svdnum,svdnum).*V(:,svdnum),'k-','Linewidth',2)
axis tight, grid on
xlabel('x_2 (km)')
set(gca,'Fontsize',20,'Linewidth',2)
subplot(122)
plot(t(plotindex)./3.15e7,S(svdnum,svdnum).*U(:,svdnum),'r-','Linewidth',4), hold on
axis tight
% xlim([0 2])
xlabel('\Deltat (yr)')

% trelax from time function (U)
funcplot = @(beta,t) beta(1).*(1-exp(-t./beta(2)));
beta0 = [0.1,1e3];
[betaval,residuals,~,covbeta] = nlinfit(t(plotindex)./3.15e7,S(svdnum,svdnum).*U(:,svdnum),funcplot,beta0);
betaci = nlparci(betaval,residuals,'covar',covbeta);
plot(t(plotindex)./3.15e7,funcplot(betaval,t(plotindex)./3.15e7),'k-')
plot(betaval(2).*[1 1],get(gca,'Ylim'),'k-','Linewidth',2)
plot(betaci(2,1).*[1 1],get(gca,'Ylim'),'k-','Linewidth',1)
plot(betaci(2,2).*[1 1],get(gca,'Ylim'),'k-','Linewidth',1)
set(gca,'Fontsize',20,'Linewidth',2)
%% estimate timeseries from 1st mode only
usvd = U(:,svdnum)*S(svdnum,svdnum)*V(:,svdnum)';
figure(12),clf
subplot(211)
imagesc(usurf(plotindex,:))
colorbar;
set(gca,'Fontsize',20,'Linewidth',2)
subplot(212)
imagesc(usvd)
colorbar;
set(gca,'Fontsize',20,'Linewidth',2)

figure(13),clf
imagesc(usurf(plotindex,:)-usvd)
colorbar;
set(gca,'Fontsize',20,'Linewidth',2)















