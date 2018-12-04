function k = getContourCurvature(x,n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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