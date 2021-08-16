function plot_imposedviscoresults(Nx,M,burger,powerval,etaval,Trecur,Vpl)

% plot results of viscocycle run
% Rishav Mallick, EOS, 2021
addpath ~/Dropbox/scripts/topotoolbox/colormaps/
addpath functions_dir/

if burger==0
    odir = ['imposedviscocycles_' num2str(Vpl,'%1.e') '/Nx_' num2str(Nx) '_M_' num2str(M) ...
    '/power_' num2str(powerval) '_' num2str(etaval) '/Trec_' num2str(Trecur)];
else
    odir = ['imposedviscocycles_' num2str(Vpl,'%1.e') '/Nx_' num2str(Nx) '_M_' num2str(M) ...
    '/burger_' num2str(etaval) '_' num2str(etaval*powerval) '/Trec_' num2str(Trecur)];
end

disp(odir)

%% load results
load([odir '/modelgeometry.mat'],'ss','shz','evl','Teq')
load([odir '/modeloutputs_f.mat'])
load([odir '/modeloutputs_s.mat'])

if isfield(shz,'L')
    x2c = shz.xx2c(:);
    x3c = shz.xx3c(:);
elseif isfield(shz,'tri')
    x2c = mean([shz.A(:,1),shz.B(:,1),shz.C(:,1)],2);
    x3c = mean([shz.A(:,2),shz.B(:,2),shz.C(:,2)],2);
end

%% plots
Teq = Trecur.*3.15e7;
ncycles = ceil(max(t)./Teq);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                    PLOTTING                          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

cspec=[cmap('steelblue',100,10,10);flipud(cmap('orange',100,57,10));flipud(cmap('orangered',100,20,25))];
% % % % % % % % % % % % % FAULT SLIP RATE % % % % % % % % % % % % 
figure(1),clf
set(gcf,'Position',[0.5 0.5 4 2].*500)
pcolor(t./3.15e7,ss.y3c./1e3,V'./Vpl), shading flat, box on
set(gca,'YDir','reverse','FontSize',15,'ColorScale','log','Color','none')
colormap(cspec)
caxis(10.^[-1 2])
ylabel('\zeta_d (km)'), xlabel('Time (yrs)')
cb=colorbar;cb.Location = 'northoutside';cb.Label.String = 'V/V_{pl}';
xlim([Teq*(ncycles-3) max(t)]./3.15e7)
print([odir '/Vtimeseries'],'-djpeg','-r200')
close

%% displacement greensfunctions
ox = linspace(0.01e3,120e3,200)';
obs = [ox,zeros(length(ox),1)];
Gd = compute_displacementkernels(obs,ss,shz);

% velocity (pin shallow creep)
V_mod = V;
V_mod(isnan(V)) = 0;

if isfield(shz,'L')
    deepx3 = max(shz.xx3c(:));
elseif isfield(shz,'tri')
    deepx3 = max(max([shz.A(:,2) shz.B(:,2) shz.C(:,2)]));
end
% deepx3 = 20e3;

% surface velocties
vsurf_f = (Gd.kd*V_mod')';
vsurf_deep = (repmat(Vpl/pi.*atan2(ox,deepx3),1,length(t)))';
vsurf_v = (Gd.l12d*e12d' +  Gd.l13d*e13d')'; % no deep loading
vsurf = 0*vsurf_f + vsurf_v + vsurf_deep;

%% rapid postseismic
figure(10),clf
Tevent = (ncycles-1)*Teq;% last cycle
tplotvec = Tevent + [1/10e3,1/5e3,1/2e3,1/1e3,1/500,1/200,1/100,1/50,1/20].*Teq;

cspec = cool(length(tplotvec));
set(gcf,'Position',[0 0 2 2].*500)
lgd = cell(length(tplotvec)+1,1);
for i = 1:length(tplotvec)
    index = find(abs(t-tplotvec(i))==min(abs(t-tplotvec(i))),1);
    plot(ox./1e3,vsurf(index,:)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
    lgd{i} = [num2str(round((t(index)-Tevent)/3.15e7*365)) ' days'];
end
plot(ox./1e3,1/pi.*atan2(ox,20e3),'k-','LineWidth',2)
lgd{length(tplotvec)+1} = '20 km plate';
legend(lgd,'Box','off')
xlabel('x_2 (km)'), ylabel('V/V_{pl}')
axis tight, grid on
ylim(10.^[-1 1.2])
set(gca,'FontSize',15,'LineWidth',2,'YScale','log')
print([odir '/Vsurf_early'],'-djpeg','-r200')
close
%% interseismic
figure(10),clf
Tevent = (ncycles-1)*Teq;% last cycle
tplotvec = Tevent + [1/50,1/20,.1,.2,.4,.5,.6,.7,.9].*Teq;

cspec = cool(length(tplotvec));
set(gcf,'Position',[0 0 2 2].*500)
lgd = cell(length(tplotvec),1);
p=[];
for i = 1:length(tplotvec)
    index = find(abs(t-tplotvec(i))==min(abs(t-tplotvec(i))),1);
    dummy = plot(ox./1e3,vsurf(index,:)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)); hold on
    p = [p;dummy];
    lgd{i} = [num2str(round((t(index)-Tevent)/3.15e7)) ' yrs'];
end
plot(ox./1e3,1/pi.*atan2(ox,20e3),'k-','LineWidth',2),text(70,0.44,[num2str(20) ' km'],'FontSize',20)
plot(ox./1e3,1/pi.*atan2(ox,40e3),'k--','LineWidth',1),text(70,0.36,[num2str(40) ' km'],'FontSize',20)
plot(ox./1e3,1/pi.*atan2(ox,60e3),'k--','LineWidth',1),text(70,0.3,[num2str(60) ' km'],'FontSize',20)
xlabel('x_2 (km)'), ylabel('V/V_{pl}')
axis tight, %grid on
legend(p,lgd,'Box','off')
ylim([0 1])
set(gca,'FontSize',15,'LineWidth',2)
print([odir '/Vsurf_late'],'-djpeg','-r200')
close
% fault transects
figure(2),clf
set(gcf,'Position',[0 0 3 1].*500)
for i = 1:length(tplotvec)
    index = find(abs(t-tplotvec(i))==min(abs(t-tplotvec(i))),1);
    plot(ss.y3c./1e3,V(index,:)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
end
axis tight,grid on, ylim([0 1.5])
xlabel('x_3 (km)'), ylabel('Normalised V')
set(gca,'FontSize',15,'Color','none','LineWidth',2)
print([odir '/Vfault_late'],'-djpeg','-r200')
% return
%% shear zone plot
figure(3),clf
set(gcf,'Position',[0 0 5 3].*500)

for i = 1:length(tplotvec)
    subplot(3,3,i)
    index = find(abs(t-tplotvec(i))==min(abs(t-tplotvec(i))),1);
    
    % plot the perturbed component
    %toplot = sqrt(e12d(index,:).^2 + e13d(index,:).^2)' - sqrt(shz.e12pl.^2 + shz.e13pl.^2);
    toplot = sqrt((e12d(index,:)-shz.e12pl').^2 + (e13d(index,:)-shz.e13pl').^2)';
    F = scatteredInterpolant(x2c,x3c,abs(toplot),'natural');
    x2g = linspace(0e3,100e3,1000)';
    x3g = linspace(min(shz.A(:,2)),max(shz.A(:,2)),200)';
    [X2g,X3g] = meshgrid(x2g,x3g);
    toplotint = F(X2g(:),X3g(:));
    
    %plotshz(shz,abs(toplot),1), shading flat, axis tight equal, box on
    imagesc(x2g./1e3,x3g./1e3,reshape(toplotint,length(x3g),length(x2g))), shading flat, hold on
    contour(x2g./1e3,x3g./1e3,reshape(toplotint,length(x3g),length(x2g)),...
        logspace(log10(Vpl*1e-6),log10(Vpl*1e-4),20),'r-','LineWidth',1)
    set(gca,'YDir','reverse','FontSize',15,'Color','none','LineWidth',2,'ColorScale','log')
    xlabel('x_2 (km)'),ylabel('x_3 (km)'),
    caxis(Vpl.*[1e-6 1e-4])
    xlim([0 1].*100), ylim([20 50])
    colormap(ttscm('oslo',20))
    cb=colorbar;cb.Label.String = 'perturbation strain rate';
    title(['\DeltaT = ' num2str(round((t(index)-Tevent)/3.15e7)) ' years']) 
end
print([odir '/shz_perturbstrainrate'],'-djpeg','-r200')

%% plot shz snapshots as velocity
tplotvec = Tevent + [.1,.2,.3,.4,.6,.9].*Teq;

figure(11),clf
set(gcf,'Position',[0 0 5 3].*500)
for i = 1:length(tplotvec)
    subplot(3,2,i)
    index = find(abs(t-tplotvec(i))==min(abs(t-tplotvec(i))),1);
    [x2g,x3g,Vg] = calculate_dispfromstrain(shz,20e3,e12d(index,:)');

    pcolor(x2g./1e3,x3g./1e3,Vg./Vpl), shading interp, hold on
    contour(x2g./1e3,x3g./1e3,Vg./Vpl,[0:0.05:0.45],'k-','LineWidth',.1)
    contour(x2g./1e3,x3g./1e3,Vg./Vpl,[0:0.1:0.4],'k-','LineWidth',2)
    contour(x2g./1e3,x3g./1e3,Vg./Vpl,[0.5:0.2:2.5],'r-','LineWidth',1)
    axis tight %equal
    xlabel('x_2 (km)'), ylabel('x_3 (km)')
    colormap(ttscm('tokyo',100))
    caxis([0 0.5])
    xlim([0 200]), 
    ylim([21 49])
    cb=colorbar;cb.Label.String = 'V/V_{pl}';
    title(['T since earthquake = ' num2str(round((t(index)-Tevent)/3.15e7)) ' years'])
    set(gca,'YDir','reverse','Fontsize',15,'LineWidth',2)
end
print([odir '/shz_integratedisp'],'-djpeg','-r200')

end