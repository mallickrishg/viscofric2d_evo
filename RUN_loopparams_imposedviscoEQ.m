% MAIN script to run viscoelastic earthquake-cycle simulations
% Brittle fault is coupled to a viscoelastic lower crust and upper mantle.
% The loading rate is entirely driven by fault slip i.e. if the fault slides steadily at the plate rate, then there will be steady shearing in the mantle. 
% Deviations from this steady-state stress state is due to elastic stress transfer from fault slip processes.
% 
% Rishav Mallick, EOS, 2021
% For imposed earthquake cycles

clear
Nx = 50;
M = 100;
Transition = 20e3;
burger = 1;

% power-law rheology
% eta = [1e18,3e18,7e18,1e19,5e19,1e20];
% power = [1:6];
% [etavec,powervec] = meshgrid(eta,power);

% burger's rheology
etaM = [1e18,5e18,1e19,5e19,1e20];
etaK = [5e17,1e18,5e18];
[etaMvec,etaKvec] = meshgrid(etaM,etaK);
powervec = etaKvec./etaMvec;
etavec = etaMvec;
index = etaKvec>etaMvec;
powervec(index) = [];
etavec(index) = [];

Trecur = 200;
Vpl = 1e-9;

Nruns = length(etavec(:));
%% run loop

parfor i = 1:Nruns
    etaval = etavec(i);
    powerval = powervec(i);
    
    RUN_MAIN_ViscoEQ_imposedcycles(Nx,M,Transition,burger,powerval,etaval,Trecur,Vpl);
end

%% plot results
if true    
    for i = 1:Nruns
        close all
        etaval = etavec(i);
        powerval = powervec(i);
        plot_imposedviscoresults(Nx,M,burger,powerval,etaval,Trecur,Vpl);
    end
end