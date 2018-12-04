function R = resolution_check(P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

diffs = @(P) P-circshift(P,1);
ratio = @(x,y) unique(max(sqrt(x.^2+y.^2)))/unique(min(sqrt(x.^2+y.^2)));

D = diffs(P);
if ratio(D(:,1),D(:,2)) > 2
    R = 0;
else
    R = 1;
end
end

