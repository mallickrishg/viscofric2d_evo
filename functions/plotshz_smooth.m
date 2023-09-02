function [x2g,x3g,toplotint] = plotshz_smooth(shz,toplot,toplot_padding)

x2c = mean([shz.A(:,1), shz.B(:,1), shz.C(:,1)],2);
x3c = mean([shz.A(:,2), shz.B(:,2), shz.C(:,2)],2);
% interpolation for contours
boundind = boundary([x2c;shz.A(:,1);shz.B(:,1);shz.C(:,1)],[x3c;shz.A(:,2);shz.B(:,2);shz.C(:,2)]);
dummy = [x2c;shz.A(:,1);shz.B(:,1);shz.C(:,1)];
x2pad = dummy(boundind);
dummy = [x3c;shz.A(:,2);shz.B(:,2);shz.C(:,2)];
x3pad = dummy(boundind);

toplotpad = ones(length(x2pad),1).*toplot_padding;

F = scatteredInterpolant([x2c;x2pad],[x3c;x3pad],abs([toplot;toplotpad]),'natural','linear');
x2g = linspace(-50e3,50e3,400)';
x3g = linspace(min(shz.A(:,2)),max(shz.A(:,2)),200)';
[X2g,X3g] = meshgrid(x2g,x3g);
toplotint = F(X2g(:),X3g(:));

% plotshz(shz,toplot,1), hold on

imagesc(x2g./1e3,x3g./1e3,reshape(toplotint,length(x3g),length(x2g)))
axis tight equal, box on
hold on
colormap(flipud(ttscm('turku',100)))


end