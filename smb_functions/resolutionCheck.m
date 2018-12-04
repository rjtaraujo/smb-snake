function R = resolutionCheck(P)
% This function returns a boolean R indicating whether the contour evolution
% may continue without a new discretization of the contour. If the maximum
% distance between two consecutive contour points is 2 or more times greater
% than the minimum, this function outputs 0, in order to trigger
% reparameterization.
%   inputs: P - the contour
%           
%   outputs: R - a boolean stating if we can continue without
%           reparameterization

diffs = @(P) P-circshift(P,1);
ratio = @(x,y) unique(max(sqrt(x.^2+y.^2)))/unique(min(sqrt(x.^2+y.^2)));

D = diffs(P);
if ratio(D(:,1),D(:,2)) > 2
    R = 0;
else
    R = 1;
end
end

