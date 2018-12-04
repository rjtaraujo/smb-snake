function B=SnakeInternalForceMatrix2D(nPoints,alpha,betas,gamma)
%
% B=SnakeInternalForceMatrix2D(nPoints,alpha,beta,gamma)
%
% inputs,
%   nPoints : The number of snake contour points
%   alpha : membrame energy  (first order)
%   betas : thin plate energy (second order)
%   gamma : Step Size (Time)
%
% outputs,
%   B : The Snake Smoothness regulation matrix
%
% Function is written by D.Kroon University of Twente (July 2010)
% and adapted by Ricardo Araújo of INESC TEC (March 2017)

% Penta diagonal matrix, one row:
alphas = ones(nPoints,1)*alpha;
betas_m1=circshift(betas,1);
betas_p1=circshift(betas,-1);

% Make the penta matrix (for every contour point)
A=repmat(betas_m1,[1, nPoints]).*circshift(eye(nPoints),2);
A=A+repmat(-alphas-2*betas_m1-2*betas,[1, nPoints]).*circshift(eye(nPoints),1);
A=A+repmat(2*alphas+betas_m1+4*betas+betas_p1,[1, nPoints]).*circshift(eye(nPoints),0);
A=A+repmat(-alphas-2*betas-2*betas_p1,[1, nPoints]).*circshift(eye(nPoints),-1);
A=A+repmat(betas_p1,[1, nPoints]).*circshift(eye(nPoints),-2);

% Calculate the inverse
B=inv(A + gamma.* eye(nPoints));


