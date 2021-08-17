function plotpatch(ss,val,kmt)

if exist('kmt','var')==1 && kmt==1
    scf=1e-3;
else
    scf = 1;
end

if isfield(ss,'Lf')
    dx = ss.Lf/2;
else
    dx = 0.1e3;
end
x = [ss.y2c-dx,ss.y2c-dx,ss.y2c+dx,ss.y2c+dx,ss.y2c-dx]';
y = [ss.y3f,ss.y3f+ss.Wf,ss.y3f+ss.Wf,ss.y3f,ss.y3f]';

patch(x.*scf,y.*scf,val)


end