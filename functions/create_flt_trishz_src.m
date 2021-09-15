function [ss,shz,src] = create_flt_trishz_src(y2i,y3i,Transition,M,Nx,scf,LS,WS,Transitionshz)

% set initial points for (top of) faults y2i,y3i
% number of patches per fault: M
% Fault depth (BDT) - Transition (m)
% shear zone discretization: Nx, scf - default values Nx = 0.2, scf = 1
% length scale(x3): LS
% length scale(x2): WS
% returns 
% ss - flt, shz - shear zone , src - deep shear zone
% structures are constructed for a fixed elastic lid thickness (Thickness ~ 20 km)
% Rishav Mallick, 2019, EOS
addpath meshing/
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                        M E S H                       %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % %  %
% Brittle-Ductile Tranisition Depth (m)
%Transition = 20e3;
ss.y3f = zeros(length(y2i)*M,1);
ss.y2f = ss.y3f;
ss.Wf = ss.y3f;
ss.M = length(y2i)*M;

for i = 1:length(y2i)
    % Fault Meshes
    
    % top of fault patches
    y3f = linspace(y3i(i),Transition,M+1)';
    ss.y3f(1+(i-1)*M:i*M) = y3f(1:end-1);
    
    ss.y2f(1+(i-1)*M:i*M) = y2i(i).*ones(M,1);
    
    % width of fault patches
    ss.Wf(1+(i-1)*M:i*M)     = diff(y3f);
    
end

% mid of fault patches
ss.y3c = ss.y3f + ss.Wf./2;
ss.y2c = ss.y2f;

%% Shear Zone Mesh

% %nx2v = 30;%scf = 20;
% dx = 0.5/Nx;%scf = 1.07;
% nx2 = round(WS/LS);
% nx3 = 2;
% dxmax = 0.05;
% % x2v = [0,logspace(0,1.6,nx2v)].*1e3;
% % x3v = [0,logspace(0,1.6,nx3v)].*1e3;
% % xvec = cumsum(tanh(linspace(0,1,nx2v)./scf));
% xvec = [];
% i = 1;
% xvec(i) = 0;
% while xvec(i)<1
%     i = i+1;
%     dxval = dx*(scf^i);
%     if dxval>dxmax
%         dxval=dxmax;
%     end
%     xvec(i) = xvec(i-1) + dxval;
% end
% xvec = (xvec-min(xvec))./(max(xvec)-min(xvec));
% 
% x2v = [xvec./3, linspace(max(xvec./3),1,round(length(xvec)/3))].*WS;
% % x3v = Transitionshz + xvec.*LS;
% % x3v = linspace(Transitionshz,Transitionshz+LS,LS./1e3);
% 
% x3dis = (1-2/pi.*atan(linspace(10,0,1.5*LS./1e3)));
% x3dis = (x3dis-min(x3dis))./(max(x3dis)-min(x3dis));
% x3v = Transitionshz + x3dis.*LS;
% 
% x2 = [x2v,                    x2v(end).*ones(1,nx3),                linspace(x2v(end),x2v(1),nx2),      x2v(1).*ones(size(x3v))]';
% x3 = [x3v(1).*ones(size(x2v)),linspace(x3v(1),x3v(end),nx3),        x3v(end).*ones(1,nx2),          fliplr(x3v)]';
% 
% p = [x2 x3];
% edge = [];
% 
p = [0,Transitionshz;...
    WS,Transitionshz;...
    WS,Transitionshz+LS;...
    0,Transitionshz+LS]./1e3;

% % for i = 1:length(x2)
% %     if i ~= length(x2)
% %         edge(i,1:2) = [i i+1];
% %     else
% %         edge(i,1:2) = [i 1];
% %     end
% % end
% % 
% % % build lowest resolution possible
% % % [vert,etri,tria,tnum] = refine2(p,edge) ;
% % 
% % hfun = 10e3 ;            % uniform "target" edge-lengths
% % % hfun = 5e3 ;            % uniform "target" edge-lengths
% % [vert,etri,tria,tnum] = refine2(p,edge,[],[],hfun) ;

x0 = 0;
y0 = Transitionshz./1e3;
hmin = Nx;
dh = scf*hmin;

hfun1 = @(x,y) hmin + dh*sqrt(( x - x0 ).^2  + ( y - y0 ).^2);

hdata = [];
hdata.fun = hfun1;

[vert,tria] = mesh2d(p, [], hdata);
close
vert = vert.*1e3;

% define triangles as A,B,C - 3 vertices (x2,x3)
shz.A = [vert(tria(:,1),1),vert(tria(:,1),2);[-vert(tria(:,1),1),vert(tria(:,1),2)]];
shz.B = [vert(tria(:,2),1),vert(tria(:,2),2);[-vert(tria(:,2),1),vert(tria(:,2),2)]];
shz.C = [vert(tria(:,3),1),vert(tria(:,3),2);[-vert(tria(:,3),1),vert(tria(:,3),2)]];

shz.tri = [tria;tria];
shz.vert = [vert;[-vert(:,1),vert(:,2)]];

shz.N = length(tria(:,1)).*2;
%% deep SOURCE
src.W = 10000e3;
src.L = 1*(2*WS);
src.xx3c = Transitionshz+LS+src.W/2;
src.xx2c = 0;

end