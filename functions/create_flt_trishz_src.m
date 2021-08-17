function [ss,shz,src] = create_flt_trishz_src(y2i,y3i,Transition,M,Nx,scf,LS,WS,Transitionshz)

% set initial points for (top of) faults y2i,y3i
% number of patches per fault: M
% Fault depth (BDT) - Transition (m)
% number of shz zones: Nx,Nz
% length scale: LS
% returns 
% ss - flt, shz - shear zone , src - deep shear zone
% structures are constructed for a fixed elastic lid thickness (Thickness ~ 20 km)
% Rishav Mallick, 2019, EOS

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

%nx2v = 30;%scf = 20;
dx = 0.5/Nx;%scf = 1.07;
nx2 = 5;nx3 = 5;

% x2v = [0,logspace(0,1.6,nx2v)].*1e3;
% x3v = [0,logspace(0,1.6,nx3v)].*1e3;
% xvec = cumsum(tanh(linspace(0,1,nx2v)./scf));
xvec = [];
i = 1;
xvec(i) = 0;
while xvec(i)<1
    i = i+1;
    xvec(i) = xvec(i-1) + dx*(scf^i);
end


x2v = xvec.*WS;
x3v = Transitionshz + xvec.*LS;


x2 = [x2v,                    x2v(end).*ones(1,nx3),                linspace(x2v(end),x2v(1),nx2),      x2v(1).*ones(size(x3v))]';
x3 = [x3v(1).*ones(size(x2v)),linspace(x3v(1),x3v(end),nx3),        x3v(end).*ones(1,nx2),          fliplr(x3v)]';

p = [x2 x3];
edge = [];
for i = 1:length(x2)
    if i ~= length(x2)
        edge(i,1:2) = [i i+1];
    else
        edge(i,1:2) = [i 1];
    end
end
% build lowest resolution possible
% [vert,etri,tria,tnum] = refine2(p,edge) ;

hfun = 10e3 ;            % uniform "target" edge-lengths
% hfun = 5e3 ;            % uniform "target" edge-lengths

[vert,etri,tria,tnum] = refine2(p,edge,[],[],hfun) ;

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