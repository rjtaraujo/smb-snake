function Bf = get_betas(Pf,Pi,Bi)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

tmp = zeros(1,2,length(Pf));
tmp(1,:,:) = Pf';

TMP = bsxfun(@(x,y) (x-y).^2, Pi, tmp);

A = squeeze(sqrt(sum(TMP,2)));

[C,I] = min(A);

Bf = Bi(I);

end

