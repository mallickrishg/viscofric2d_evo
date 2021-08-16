function u1=computeDisplacementAntiplaneTriangleShearZone( ...
    x2,x3,A,B,C,e12,e13)
% function COMPUTEDISPLACEMENTANTIPLANETRIANGLESHEARZONE computes the
% displacement field associated with deforming triangular shear zones
% considering the following geometry:
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
% e12, e13           source strain component 12 and 13 in the shear zone.
%
% Output:
% u1                 displacement component in the north (along-strike)
%                    direction.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - March 4, 2018, Singapore.

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

u1=2*(e12*nC(1)+e13*nC(2))*T11(A,B) ...
  +2*(e12*nA(1)+e13*nA(2))*T11(B,C) ...
  +2*(e12*nB(1)+e13*nB(2))*T11(C,A);

    function u=T11(A,B)
        % radial vector
        a=(B-A)/2;
        a=a/norm(a);
        % normal vector
        n=[-a(2);a(1)];
        
        Gamma=@(r) -(a(1)*(x2-r(1))+a(2)*(x3-r(2))) .* log( (x2-r(1)).^2+(x3-r(2)).^2 ) ...
                 -2*(n(1)*(x2-r(1))+n(2)*(x3-r(2))) .* atan( (a(1)*(x2-r(1))+a(2)*(x3-r(2))) ./ ...
                                                             (n(1)*(x2-r(1))+n(2)*(x3-r(2))) ) ;
                                         
        u=-(Gamma(B)-Gamma(A)) / (4*pi);
        
        % image
        A(2)=-A(2);
        B(2)=-B(2);
        
        % radial vector
        a=(B-A)/2;
        a=a/norm(a);
        % normal vector
        n=[-a(2);a(1)];
        
        Gamma=@(r) -(a(1)*(x2-r(1))+a(2)*(x3-r(2))) .* log( (x2-r(1)).^2+(x3-r(2)).^2 ) ...
                 -2*(n(1)*(x2-r(1))+n(2)*(x3-r(2))) .* atan( (a(1)*(x2-r(1))+a(2)*(x3-r(2))) ./ ...
                                                             (n(1)*(x2-r(1))+n(2)*(x3-r(2))) ) ;
                                         
        u=u-(Gamma(B)-Gamma(A)) / (4*pi);
    end

end











