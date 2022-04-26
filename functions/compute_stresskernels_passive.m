function evl = compute_stresskernels_passive(obs,src)

Nobs = length(obs(:,1));

G = 30e3; %MPa

% Stress kernels for fault slip
s12h=@(x2,x3,y2,y3,Wf) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

s13h=@(x2,x3,y2,y3,Wf) G*( ...
    (x2-y2)./((x2-y2).^2+(x3-y3).^2)-(x2-y2)./((x2-y2).^2+(x3+y3).^2) ...
    -(x2-y2)./((x2-y2).^2+(x3-y3-Wf).^2)+(x2-y2)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

% Stress kernels from fault
if isfield(src,'y2c')
    ss = src;
    evl.k12=zeros(Nobs,ss.M);
    evl.k13=zeros(Nobs,ss.M);
    
    for k=1:ss.M
        evl.k12(:,k)=s12h(obs(:,1),obs(:,2),ss.y2f(k),ss.y3f(k),ss.Wf(k));
        evl.k13(:,k)=s13h(obs(:,1),obs(:,2),ss.y2f(k),ss.y3f(k),ss.Wf(k));
        
    end
elseif isfield(src,'tri')
    shz = src;
    evl.l1212=zeros(Nobs,shz.N);
    evl.l1213=zeros(Nobs,shz.N);
    evl.l1312=zeros(Nobs,shz.N);
    evl.l1313=zeros(Nobs,shz.N);
    
    for k=1:(shz.N)
        [s12,s13] = computeStressAntiplaneTriangleShearZone(obs(:,1),obs(:,2),shz.A(k,:),shz.B(k,:),shz.C(k,:),1,0,G);
        evl.l1212(:,k) = s12;
        evl.l1213(:,k) = s13;
        
        [s12,s13] = computeStressAntiplaneTriangleShearZone(obs(:,1),obs(:,2),shz.A(k,:),shz.B(k,:),shz.C(k,:),0,1,G);
        evl.l1312(:,k) = s12;
        evl.l1313(:,k) = s13;
    end
else
    error('invalid source object')
end




end