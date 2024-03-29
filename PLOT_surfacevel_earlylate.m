% plot early (peak) postseismic and late interseismic velocities with a linear scale
% Rishav Mallick, EOS, 2021

clear
addpath functions/

% specify parameters
Nx = 0.2;
M = 100;
Vpl = 1e-9;
Trecur = 50;

% specify viscosity/rheological param and power exponent
etavec = [3e18,5e19,3e18];
powervec = [1,5e18/etavec(2),3];
burgervec = [0,1,0];


% specify surface observation points
ox = linspace(-100e3,100e3,400)';

%% plot results
figure(1),clf
set(gcf,'Position',[0 0.5 1.5 1.5].*500)
% cspec = [rgb('orangered');rgb('sky blue');rgb('brown')];
cspec = [1,0,0;0,1,0;0,0,1];
for count = 1:length(etavec)
    
    etaval = etavec(count);
    powerval = powervec(count);
    burger = burgervec(count);
    
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
    
    
    % PLOT surface velocities at specified times since the earthquake
    tplotvec = [0.01*3.15e7;Teq];
    
    % compute surface velocities    
    obs = [ox,zeros(length(ox),1)];
    Gd = compute_displacementkernels(obs,ss,shz);
    deepx3 = max(shz.A(:,2));
    vsurf_deep = (repmat(Vpl/pi.*atan2(ox,deepx3),1,length(t)))';
    vsurf_v = (Gd.l12d*e12d' +  Gd.l13d*e13d')'; % no deep loading
    vsurf = vsurf_v + vsurf_deep;
    
    for i = 1:length(tplotvec)
        index = find(abs(t-tplotvec(i))==min(abs(t-tplotvec(i))),1);        
        plot(ox./1e3,vsurf(index,:)./Vpl,'-','Linewidth',1.5,'Color',cspec(count,:)), hold on
    end

end
plot(ox./1e3,1/pi*atan2(ox,20e3),'k--','Linewidth',1)
axis tight
ylabel('v/v^{\infty}'),xlabel('x_2 (km)')
legend('Linear Maxwell','','Burgers','','Power-law','','Elastic plate')
set(legend,'box','off')
ylim([-1 1]*1.8)
set(gca,'FontSize',20,'LineWidth',1,'YScale','lin')
% xlim([-1 1]*200)

%% plot as 3 separate plots with various snapshots
figure(10),clf
set(gcf,'Position',[0 0 4 1.1]*500,'Color','w')
for count = 1:length(etavec)
    %figure(count+1),clf
    %set(gcf,'Position',[(count-1)*1.5 0 1.5 1.5].*500)
    etaval = etavec(count);
    powerval = powervec(count);
    burger = burgervec(count);
    
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
    
    
    
    % PLOT surface velocities at specified times
    %tplotearly = logspace(-3,1,10)'.*3.15e7;
    %tplotlate = [0.1,0.3,0.5,0.7,0.9]'*Teq;
    
    %tplotvec = [tplotearly;tplotlate];
    tplotvec = [0.001,0.01,0.1,0.5,0.9]'*Teq;
    
    % compute surface velocities    
    Gd = compute_displacementkernels(obs,ss,shz);
    deepx3 = max(shz.A(:,2));
    vsurf_deep = (repmat(Vpl/pi.*atan2(ox,deepx3),1,length(t)))';
    vsurf_v = (Gd.l12d*e12d' +  Gd.l13d*e13d')'; % no deep loading
    vsurf = vsurf_v + vsurf_deep;
    
    figure(10)
    subplot(1,length(etavec),count)
    %cspec = [cool(length(tplotearly));parula(length(tplotlate))];
    cspec = cool(length(tplotvec));
    for i = 1:length(tplotvec)
        index = find(abs(t-tplotvec(i))==min(abs(t-tplotvec(i))),1);
        
        plot(ox./1e3,vsurf(index,:)./Vpl,'-','Linewidth',3,'Color',cspec(i,:)), hold on
    end
    plot(ox./1e3,1/pi*atan2(ox,20e3),'k:','Linewidth',2)
    axis tight
    ylim([-1 1]*2)
    ylabel('v/v^{\infty}')
    
    if count==1
        title('Linear Maxwell','FontWeight','normal')
    elseif count==2
        title('Linear Burgers','FontWeight','normal')
        legend('\Deltat/T_{eq} = 0.001','\Deltat/T_{eq} = 0.01','\Deltat/T_{eq} = 0.1','\Deltat/T_{eq} = 0.5','\Deltat/T_{eq} = 0.9','steady interseismic',...
        'location','southeast','box','off')
    else
        title('Power Law','FontWeight','normal')
    end    
    
    set(gca,'FontSize',20,'LineWidth',1,'YScale','lin','YTick',[-2:1:2])
    xlim([-1 1]*100)    
    xlabel('x_2 (km)')

    figure(20+count),clf
    set(gcf,'Position',[0 0 3 1.4]*500,'Color','w')
    for i = 1:length(tplotvec)
        index = find(abs(t-tplotvec(i))==min(abs(t-tplotvec(i))),1);        
        plot(ox./1e3,vsurf(index,:)./Vpl,'-','Linewidth',4,'Color',cspec(i,:)), hold on
    end
    plot(ox./1e3,1/pi*atan2(ox,20e3),'k-','Linewidth',4)
    axis tight
    ylim([-1 1]*2)
    ylabel('v/v^{\infty}')
    set(gca,'FontSize',40,'LineWidth',2,'YScale','lin','YTick',[-2:1:2],'TickDir','out')
    xlim([-1 1]*100)    
    xlabel('x_2 (km)')
    legend('\Deltat/T_{eq} = 0.001','\Deltat/T_{eq} = 0.01','\Deltat/T_{eq} = 0.1','\Deltat/T_{eq} = 0.5','\Deltat/T_{eq} = 0.9',...
        ...
        'location','northwest','box','off')
    print(['Figures/postseismic_comparerheologies_' num2str(count)],'-djpeg','-r300')

end
%print(['Figures/postinterseismic_comparerheologies_'],'-djpeg','-r300')
