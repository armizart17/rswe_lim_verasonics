
% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_TV(b,A,mu,M,N,tol,mask,minimask)
% function u = IRLS_TV(b,A,mu,M,N,tol,mask,minimask)

AtA = A'*A;
Atb = A'*b;

%figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dx = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

%[u,~] = cgs(AtA+mu*(D')*D,Atb);
[u,~] = minres(AtA,Atb);
G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask,1,false);

while error > tol && ite_irls < 200
    
    X = reshape(u,M,N);
    ite_irls = ite_irls + 1;
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];
    
    vksquare = Dh.^2 + Dv.^2;
    vksquare = vksquare(:);
    
    % eps = 0.3;
    eps = 0.3;
    P = sqrt(vksquare + eps^2);
    P = 1./P;
    %Huber seems not to work;
    
    %P = sqrt(P.^2 + eps^2);
    P = P(:).*minimask;
    omega = speye(M*N);   % sparse diagonal of just ones
    omega = spdiags(P,0,omega);   % sparse diagonal of P values instead of ones.
    W = kron(speye(2),omega);
    
    [u,~] = minres(AtA + mu*D'*W*D, Atb, tol, 200);
    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask, 1, false);
%     G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u,M,N,minimask);
    % error = abs(G(ite_irls+1) - G(ite_irls));
    error = (G(ite_irls+1) - G(ite_irls))/G(ite_irls);
    
end
% figure,plot(G)
%figure(909); plot(1:length(G),G);

end

function [TV] = TVcalc(B,M,N,mask,isotropic,colMajor)
% Returns the Total Variation
% Inputs: 
%       B               vector containing an image
%       M,N             image size
%       mask            binary mask of analyzed region
%       isotropic       true if the TV is isotropic, false if anisotropic
%       colMajor        true if the matrix is stored col-major, false if
%                       stored by rows
%  
% Outputs:
%       TV              number
%
% Author: Andres Leonel Coila
% Modified by Sebastian Merino

mask(isnan(mask)) = 0;
mask = mask(:);

if colMajor
    X = reshape(B,M,N);
else
    X = reshape(B,N,M)';
end

Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

if (isotropic)
    P = Dh.^2 + Dv.^2;
    P = sqrt(P);
    TV = norm(P(:).*mask,1);
else
    P = abs(Dh) + abs(Dv);
    TV = norm(P(:).*mask,1);
end

end