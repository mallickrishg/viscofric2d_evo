clear
% close all

x0 = 0;
y0 = 20;
hfun1 = @(x,y) 0.1 + 0.3*sqrt(( x - x0 ).^2  + ( y - y0 ).^2);

% v =[0.0, 0.0;...
%     2.0, 0.0;...
%     2.0, 1.0;...
%     2.0, 2.0;...
%     1.0, 2.0;...
%     0.0, 2.0 ];

v = [0,20;...
    210,20;...
    210,50;...
    0,50];

hdata = [];
hdata.fun = hfun1;


[ p, t ] = mesh2d ( v, [], hdata);
close
figure(1),clf
plot(v(:,1),v(:,2),'ko','MarkerFaceColor','r')
hold on
for i = 1:length(t)
    plot([p(t(i,:),1);p(t(i,1),1)],[p(t(i,:),2);p(t(i,1),2)],'r-')
end
axis tight equal
set(gca,'YDir','reverse','FOntsize',20,'Linewidth',2)