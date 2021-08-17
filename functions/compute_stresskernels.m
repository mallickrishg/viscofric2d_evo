function evl = compute_stresskernels(ss,shz,src)
% compute stress kernels for flt, shz using analytical expressions derived
% in Lambert & Barbot, 2016 GRL
% returns evl structure with all apropriate kernels
% Rishav Mallick, EOS, 2019

G = 30e3; %MPa

% Boxcar function
boxc=@(x) (x+0.5>=0)-(x-0.5>=0);

% Stress kernels for distributed deformation (source(i)-receiver(j))
s1312=@(D,L,W,x2,x3) G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    +log((x2-L/2).^2+(x3+D+W).^2) - log((x2+L/2).^2+(x3+D+W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2) ...
    -log((x2-L/2).^2+(x3+D).^2) + log((x2+L/2).^2+(x3+D).^2));

s1212=@(D,L,W,x2,x3) G/pi*( ...
    atan((x3-D)./(x2+L/2))-atan((x3-D)./(x2-L/2)) ...
    +atan((x3-D-W)./(x2-L/2))-atan((x3-D-W)./(x2+L/2)) ...
    -atan((x3+D+W)./(x2-L/2))-atan((x3+D)./(x2+L/2)) ...
    +atan((x3+D)./(x2-L/2))+atan((x3+D+W)./(x2+L/2)))...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);

s1213=@(D,L,W,x2,x3) G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    -log((x2-L/2).^2+(x3+D+W).^2) + log((x2+L/2).^2+(x3+D+W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2) ...
    +log((x2-L/2).^2+(x3+D).^2) - log((x2+L/2).^2+(x3+D).^2));

s1313=@(D,L,W,x2,x3) G/pi*( ...
    atan((x2+L/2)./(x3-D))  -atan((x2-L/2)./(x3-D)) ...
    -atan((x2+L/2)./(x3-D-W))+atan((x2-L/2)./(x3-D-W)) ...
    +atan((x2+L/2)./(x3+D))  -atan((x2-L/2)./(x3+D)) ...
    -atan((x2+L/2)./(x3+D+W))+atan((x2-L/2)./(x3+D+W)))...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);

% Stress kernels for fault slip
s12h=@(x2,x3,y2,y3,Wf) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

s13h=@(x2,x3,y2,y3,Wf) G*( ...
    (x2-y2)./((x2-y2).^2+(x3-y3).^2)-(x2-y2)./((x2-y2).^2+(x3+y3).^2) ...
    -(x2-y2)./((x2-y2).^2+(x3-y3-Wf).^2)+(x2-y2)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;


%% Create stress kernels for fault and shear zone interactions

% Stress kernels from fault
evl.kl12=zeros(shz.N,ss.M);
evl.kl13=zeros(shz.N,ss.M);

evl.K=zeros(ss.M,ss.M);   % Fault self stress

% Fields from Faults
% Evaluate the stress at the center of the fault and shear zone patches
% and displacement at the top of the patch

if isfield(ss,'Lf')
    for k=1:ss.M
        % Stress on shear zones from fault slip
        evl.kl12(:,k)=s1212(ss.y3f(k),ss.Lf(k),ss.Wf(k),shz.xx2c(:)-ss.y2c(k),shz.xx3c(:));
        evl.kl13(:,k)=s1213(ss.y3f(k),ss.Lf(k),ss.Wf(k),shz.xx2c(:)-ss.y2c(k),shz.xx3c(:));
        
        % Stress on faults from fault slip
        evl.K(:,k) =  s1212(ss.y3f(k),ss.Lf(k),ss.Wf(k),ss.y2c(:)-ss.y2c(k),ss.y3c(:));
    end
    
else    
    for k=1:ss.M
        % Stress on shear zones from fault slip
        if isfield(shz,'L')
            evl.kl12(:,k)=s12h(shz.xx2c(:),shz.xx3c(:),ss.y2f(k),ss.y3f(k),ss.Wf(k));
            evl.kl13(:,k)=s13h(shz.xx2c(:),shz.xx3c(:),ss.y2f(k),ss.y3f(k),ss.Wf(k));
        elseif isfield(shz,'tri')
            x2c = mean([shz.A(:,1),shz.B(:,1),shz.C(:,1)],2);
            x3c = mean([shz.A(:,2),shz.B(:,2),shz.C(:,2)],2);
            evl.kl12(:,k)=s12h(x2c(:),x3c(:),ss.y2f(k),ss.y3f(k),ss.Wf(k));
            evl.kl13(:,k)=s13h(x2c(:),x3c(:),ss.y2f(k),ss.y3f(k),ss.Wf(k));
        else
            disp('Not a valid shz structure')
        end
        
        % Stress on faults from fault slip
        evl.K(:,k)=s12h(ss.y2f,ss.y3c,ss.y2f(k),ss.y3f(k),ss.Wf(k));
    end
end

% Stress kernels from shear zones
evl.l1312=zeros(shz.N); 
evl.l1313=zeros(shz.N);
evl.l1212=zeros(shz.N);
evl.l1213=zeros(shz.N);

evl.lk1212=zeros(length(ss.y3f),shz.N);
evl.lk1312=zeros(length(ss.y3f),shz.N);

% Fields from Shear zones
% Evaluate stress at the center of the fault and shear zone patches
% and displacement at the top of the patch
for k=1:(shz.N)
    % Stress on shear zones due to strain in shear zones
    if isfield(shz,'L')
        evl.l1212(:,k) = s1212(shz.xx3c(k)-shz.W(k)/2,shz.L(k),shz.W(k),shz.xx2c(:)-shz.xx2c(k),shz.xx3c(:));
        evl.l1312(:,k) = s1312(shz.xx3c(k)-shz.W(k)/2,shz.L(k),shz.W(k),shz.xx2c(:)-shz.xx2c(k),shz.xx3c(:));
        evl.l1213(:,k) = s1213(shz.xx3c(k)-shz.W(k)/2,shz.L(k),shz.W(k),shz.xx2c(:)-shz.xx2c(k),shz.xx3c(:));
        evl.l1313(:,k) = s1313(shz.xx3c(k)-shz.W(k)/2,shz.L(k),shz.W(k),shz.xx2c(:)-shz.xx2c(k),shz.xx3c(:));
        
        % Stress on faults due to strain in shear zones
        evl.lk1212(:,k) = s1212(shz.xx3c(k)-shz.W(k)/2,shz.L(k),shz.W(k),ss.y2c(:)-shz.xx2c(k),ss.y3c(:));
        evl.lk1312(:,k) = s1312(shz.xx3c(k)-shz.W(k)/2,shz.L(k),shz.W(k),ss.y2c(:)-shz.xx2c(k),ss.y3c(:));
    elseif isfield(shz,'tri')
        %x2c = mean([shz.A(:,1),shz.B(:,1),shz.C(:,1)],2);
        %x3c = mean([shz.A(:,2),shz.B(:,2),shz.C(:,2)],2);
        % [s12,s13] = computeStressAntiplaneTriangleShearZone(x2,x3,A,B,C,e12,e13,G);
        
        [s12,s13] = computeStressAntiplaneTriangleShearZone(x2c,x3c,shz.A(k,:),shz.B(k,:),shz.C(k,:),1,0,G);
        evl.l1212(:,k) = s12;
        evl.l1213(:,k) = s13;
        
        [s12,s13] = computeStressAntiplaneTriangleShearZone(x2c,x3c,shz.A(k,:),shz.B(k,:),shz.C(k,:),0,1,G);
        evl.l1312(:,k) = s12;
        evl.l1313(:,k) = s13;
        
        % Stress on faults due to strain in shear zones
        [s12,~] = computeStressAntiplaneTriangleShearZone(ss.y2c,ss.y3c,shz.A(k,:),shz.B(k,:),shz.C(k,:),1,0,G);
        evl.lk1212(:,k) = s12;
        [s12,~] = computeStressAntiplaneTriangleShearZone(ss.y2c,ss.y3c,shz.A(k,:),shz.B(k,:),shz.C(k,:),0,1,G);
        evl.lk1312(:,k) = s12;
        
    else
        disp('Not a valid shz structure')
    end
    
end

%% If src exists
if nargin==3
    % Stress Field from src to fault
    evl.srck1212 = s1212(src.xx3c-src.W/2,src.L,src.W,ss.y2c(:)-src.xx2c,ss.y3c(:));
    evl.srck1312 = s1312(src.xx3c-src.W/2,src.L,src.W,ss.y2c(:)-src.xx2c,ss.y3c(:));
    
    % Stress Field from src to shz
    evl.srcl1212 = s1212(src.xx3c-src.W/2,src.L,src.W,shz.xx2c(:)-src.xx2c,shz.xx3c(:));
    evl.srcl1312 = s1312(src.xx3c-src.W/2,src.L,src.W,shz.xx2c(:)-src.xx2c,shz.xx3c(:));
    evl.srcl1213 = s1213(src.xx3c-src.W/2,src.L,src.W,shz.xx2c(:)-src.xx2c,shz.xx3c(:));
    evl.srcl1313 = s1313(src.xx3c-src.W/2,src.L,src.W,shz.xx2c(:)-src.xx2c,shz.xx3c(:));
end




end