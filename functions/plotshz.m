function plotshz(shz,val,kmt)


if exist('kmt','var')==1 && kmt==1
    scf = 1/1e3;
else 
    scf = 1;
end

if isfield(shz,'L')
    dx = shz.L(:)/2;
    dy = shz.W(:)/2;
    
    x = [shz.xx2c(:)-dx,shz.xx2c(:)-dx,shz.xx2c(:)+dx,shz.xx2c(:)+dx,shz.xx2c(:)-dx]';
    y = [shz.xx3c(:)-dy,shz.xx3c(:)+dy,shz.xx3c(:)+dy,shz.xx3c(:)-dy,shz.xx3c(:)-dy]';
    
    patch(x.*scf,y.*scf,val)

elseif isfield(shz,'tri')
    x = [shz.A(:,1),shz.B(:,1),shz.C(:,1)]';
    y = [shz.A(:,2),shz.B(:,2),shz.C(:,2)]';
    
    patch(x.*scf,y.*scf,val)
    
else
    disp('Nothing to plot')
end

end