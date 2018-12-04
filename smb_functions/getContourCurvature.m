function k = getContourCurvature(x,n)
% This function returns a vector k comprising the local curvature at each
% point of the contour x. The curvature at point x_i is calculated as the 
% inverse of the radius of the circle that better fits x_(i-n):x_(i+n)
%   inputs: x - the contour
%           n - the number of neighbours taken into account in the circle
%               fitting process
%           
%   outputs: k - curvature along the contour

k = zeros(1,length(x));
x_row = x(:,1);
x_col = x(:,2);

lastwarn('');

for p = 1 : size(x,1)

    N = -n:n;
    
    x_row_p = x_row(1+mod(p-1-N,length(x)));
    x_col_p = x_col(1+mod(p-1-N,length(x)));
    
	m_row = mean(x_row_p); m_col = mean(x_col_p);
    X = x_row_p - m_row; Y = x_col_p - m_col; % Get differences from means
    
    dx2 = mean(X.^2); dy2 = mean(Y.^2); % Get variances
    t = [X,Y]\(X.^2-dx2+Y.^2-dy2)/2; % Solve least mean squares problem
    
    [~,msgid] = lastwarn;
    if isequal(msgid,'MATLAB:rankDeficientMatrix')
        k(p) = 0;
        lastwarn('');
        continue
    end
    a0 = t(1); b0 = t(2); % t is the 2 x 1 solution array [a0;b0]
    r = sqrt(dx2+dy2+a0^2+b0^2); % Calculate the radius
    k(p) = 1/r; % Get the curvature
end