function Bf = getBetas(Pf,Pi,Bi)
% This function associates the rigidity coefficients Bi to the newly
% reparameterized contour Pf.
%   inputs: Pf - the reparameterized contour
%           Pi - the contour before reparameterization
%           Bi - rigidity coefficients of Pi
%           
%   outputs: Bf - rigidity coefficients of Pf

tmp = zeros(1,2,length(Pf));
tmp(1,:,:) = Pf';

TMP = bsxfun(@(x,y) (x-y).^2, Pi, tmp);

A = squeeze(sqrt(sum(TMP,2)));

[C,I] = min(A);

Bf = Bi(I);

end

