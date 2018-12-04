function Pf = reparameterize(Px,nPoints)
% This function returns a reparameterization Pf of a contour Px,
% parameterized by nPoints
%   inputs: Px - the contour
%           nPoints - the number of discrete points of the contour
%           
%   outputs: Pf - the reparameterized contour
    
    D = sqrt(sum((circshift(Px,-1)-Px).^2,2));
    idx = [0; cumsum(D)];
    perimeter = idx(end);
    
    signal_x = [Px(:,1); Px(1,1)];
    signal_y = [Px(:,2); Px(1,2)];
    
    xq = 0:perimeter/nPoints:perimeter-perimeter/nPoints;
    Pf(:,1) = interp1(idx,signal_x,xq);
    Pf(:,2) = interp1(idx,signal_y,xq);
    
end

