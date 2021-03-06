function x = getInitialContour(c, row_rad, col_rad, n_points)
% This function returns an initial discrete ellipsoid contour given its center
% coordinates, vertical and horizontal radius length, and the number of discrete
% points to use
%   inputs: c - 2-element vector, where c(1) and c(2) are the row and
%               column of the circle's center
%           row_rad - vertical radius length
%           col_rad - horizontal radius length
%           n_points - the number of points discretizing the contour
%
%   outputs: x - the coordinates of the discrete contour

    angle_step = (2*pi)/n_points;
    angles = angle_step:angle_step:2*pi;
    ptx_row = c(1) + row_rad * sin(angles);
    ptx_col = c(2) + col_rad * cos(angles);
    x = [ptx_row(:), ptx_col(:)];

end

