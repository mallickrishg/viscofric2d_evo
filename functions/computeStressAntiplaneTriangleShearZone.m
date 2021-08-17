function [s12,s13]=computeStressAntiplaneTriangleShearZone( ...
    x2,x3,A,B,C,e12,e13,G)
% function COMPUTESTRESSANTIPLANETRIANGLESHEARZONE computes the
% stress field associated with deforming triangular shear zones
% considering the following geometry
%
%              surface
%      -------------+-------------- E (x2)
%                   |
%                   |     + A
%                   |    /  .
%                   |   /     .
%                   |  /        .
%                   | /           .
%                   |/              + B
%                   /            .
%                  /|          /
%                 / :       .
%                /  |    /
%               /   : .
%              /   /|
%             / .   :
%            +      |
%          C        :
%                   |
%                   D (x3)
%
% Input:
% x2, x3             east coordinates and depth of the observation point,
% A, B, C            east and depth coordinates of the vertices,
% e12, e13           source strain component 12 and 13 in the shear zone,
% G                  rigidity in the half space.
%
% Output:
% s12                horizontal stress component,
% s13                vertical stress component.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - March 5, 2018, Singapore.

assert(min(x3(:))>=0,'depth must be positive.');

% unit vectors
nA = [C(2)-B(2);
    B(1)-C(1)]/norm(C-B);
nB = [C(2)-A(2);
    A(1)-C(1)]/norm(C-A);
nC = [B(2)-A(2);
    A(1)-B(1)]/norm(B-A);

% check that unit vectors are pointing outward
if (nA'*(A(:)-(B(:)+C(:))/2))>0
    nA=-nA;
end
if (nB'*(B(:)-(A(:)+C(:))/2))>0
    nB=-nB;
end
if (nC'*(C(:)-(A(:)+B(:))/2))>0
    nC=-nC;
end

s12=2*(e12*nC(1)+e13*nC(2))*T112(A,B) ...
   +2*(e12*nA(1)+e13*nA(2))*T112(B,C) ...
   +2*(e12*nB(1)+e13*nB(2))*T112(C,A);

s13=2*(e12*nC(1)+e13*nC(2))*T113(A,B) ...
   +2*(e12*nA(1)+e13*nA(2))*T113(B,C) ...
   +2*(e12*nB(1)+e13*nB(2))*T113(C,A);

% remove anelastic strain
Omega=@(x2,x3) heaviside(((A(1)+B(1))/2-x2)*nC(1)+((A(2)+B(2))/2-x3)*nC(2)) ...
             .*heaviside(((B(1)+C(1))/2-x2)*nA(1)+((B(2)+C(2))/2-x3)*nA(2)) ...
             .*heaviside(((C(1)+A(1))/2-x2)*nB(1)+((C(2)+A(2))/2-x3)*nB(2));

s12=2*G*(s12-Omega(x2,x3)*e12);
s13=2*G*(s13-Omega(x2,x3)*e13);

    function y=heaviside(x)
        y=x>=0;
    end

    function u1d2=T112(A,B)
        % radial vector
        a=(B-A)/2;
        a=a/norm(a);
        % normal vector
        n=[-a(2);a(1)];
        
        Gamma2=@(r) a(1)*log( (x2-r(1)).^2+(x3-r(2)).^2 )/(8*pi) ...
                   +n(1)*atan( (a(1)*(x2-r(1))+a(2)*(x3-r(2))) ./ ...
                               (n(1)*(x2-r(1))+n(2)*(x3-r(2))) )/(4*pi);
                                                         
        u1d2=Gamma2(B)-Gamma2(A);
        
        % image
        A(2)=-A(2);
        B(2)=-B(2);
        
        % radial vector
        a=(B-A)/2;
        a=a/norm(a);
        % normal vector
        n=[-a(2);a(1)];
        
        Gamma2=@(r) a(1)*log( (x2-r(1)).^2+(x3-r(2)).^2 )/(8*pi) ...
                   +n(1)*atan( (a(1)*(x2-r(1))+a(2)*(x3-r(2))) ./ ...
                               (n(1)*(x2-r(1))+n(2)*(x3-r(2))) )/(4*pi);
        
        u1d2=u1d2+(Gamma2(B)-Gamma2(A));
    end

    function u=T113(A,B)
        % radial vector
        a=(B-A)/2;
        a=a/norm(a);
        % normal vector
        n=[-a(2);a(1)];
        
        Gamma3=@(r) a(2)*log( (x2-r(1)).^2+(x3-r(2)).^2 )/(8*pi) ...
                   +n(2)*atan( (a(1)*(x2-r(1))+a(2)*(x3-r(2))) ./ ...
                               (n(1)*(x2-r(1))+n(2)*(x3-r(2))) )/(4*pi);
        
        u=Gamma3(B)-Gamma3(A);
        
        % image
        A(2)=-A(2);
        B(2)=-B(2);
        
        % radial vector
        a=(B-A)/2;
        a=a/norm(a);
        % normal vector
        n=[-a(2);a(1)];
        
        Gamma3=@(r) a(2)*log( (x2-r(1)).^2+(x3-r(2)).^2 )/(8*pi) ...
                   +n(2)*atan( (a(1)*(x2-r(1))+a(2)*(x3-r(2))) ./ ...
                               (n(1)*(x2-r(1))+n(2)*(x3-r(2))) )/(4*pi);
        
        u=u+(Gamma3(B)-Gamma3(A));
    end

end












