function [X2,J] = SMB_Snake2D(I,X,Options)
% This is a function extending the Snake_2D function of Dirk-Jan Kroon to
% incorporate the capabilities of the SMB snake

% Process inputs
defaultoptions=struct('Verbose',false,'nPoints',100,'Wline',0.04,'Wedge',2,'Wterm',0.01,'Sigma1',10,'Sigma2',20,'Alpha',0.2,'Beta',0.2,'Delta',0.1,'Gamma',1,'Kappa',2,'Iterations',100,'GIterations',0,'Mu',0.2,'Sigma3',1,'betaMean',5,'maxValue',10,'lambda',1);
if(~exist('Options','var')), 
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))), 
        warning('snake:unknownoption','unknown options found');
    end
end

% The contour must always be clockwise (because of the balloon force)
X=MakeContourClockwise2D(X);

% Make an uniform sampled contour description
X=InterpolateContourPoints2D(X,Options.nPoints);

% Transform the Image into an External Energy Image
Eext = ExternalForceImage2D(I,Options.Wline, Options.Wedge, Options.Wterm, Options.Sigma1);

% Make the external force (flow) field.
Fx=ImageDerivatives2D(Eext,Options.Sigma2,'x');
Fy=ImageDerivatives2D(Eext,Options.Sigma2,'y');
Fext(:,:,1)=-Fx*2*Options.Sigma2^2;
Fext(:,:,2)=-Fy*2*Options.Sigma2^2;

% Do Gradient vector flow, optimalization
Fext=GVFOptimizeImageForces2D(Fext, Options.Mu, Options.GIterations, Options.Sigma3);

% Make the interal force matrix, which constrains the moving points to a
% smooth contour
S=SnakeInternalForceMatrix2D(Options.nPoints,Options.Alpha,Options.Beta,Options.Gamma);
h=[];

% Show the image, contour and force field
if(Options.Verbose)
    h = figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6]); set(h,'render','opengl')
    subplot(2,2,1),
    imshow(I,[]); hold on; plot(X(:,2),X(:,1),'b.'); hold on;
    title('Image and initial contour')
    subplot(2,2,2),
    imshow(Eext,[]); 
    title('The external energy');
    subplot(2,2,3), 
    [x,y]=ndgrid(1:2:size(Fext,1),1:2:size(Fext,2));
    imshow(I), hold on; quiver(y,x,Fext(1:2:end,1:2:end,2),Fext(1:2:end,1:2:end,1));
    title('The external force field ')
    subplot(2,2,4), 
    imshow(Eext,[]), hold on; plot(X(:,2),X(:,1),'b.'); hc1 = plot(X(:,2),X(:,1),'r.');
    title('Snake movement ')
    
    drawnow
end

for i=1:Options.Iterations    
    
    X=SnakeMoveIteration2D(S,X,Fext,Options.Gamma,Options.Kappa,Options.Delta);

    if ~resolution_check(X)
        X = reparameterize(X,Options.nPoints);
    end
    
    % Show current contour
    if(Options.Verbose)
        set(hc1, 'XData', X(:,2)', 'YData', X(:,1)');
        drawnow
    end
end

%%
% after fitting the traditional snake model, we analyze the curvature along
% the contour x, optimize the beta distribution and further evolve the
% contour, to obtain the final result x2

k = getContourCurvature(X,1);

init_betas = ones(1, Options.nPoints) * Options.betaMean;
opt_betas = optimizer_l0_m(Options.lambda, Options.maxValue, 100, k, init_betas);
opt_betas = opt_betas';

X2 = X;
betas = opt_betas; 

% Make the interal force matrix, which constrains the moving points to a
% smooth contour
S=SnakeInternalForceMatrix2D_pos(Options.nPoints,Options.Alpha,opt_betas,Options.Gamma);
h2=[];

if(Options.Verbose)
    h2 = figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6]); set(h2,'render','opengl')
    subplot(1,3,1), plot(k); title('curvature');
    subplot(1,3,2), plot(init_betas,'b'); hold on; plot(opt_betas,'r'); title('beta optimization');
    subplot(1,3,3), 
    imshow(Eext,[]), hold on; plot(X2(:,2),X2(:,1),'b.'); hc2 = plot(X2(:,2),X2(:,1),'r.');
    title('Final evolution')
    
    drawnow
    pause(1)
end

for i=1:Options.Iterations
    X2=SnakeMoveIteration2D(S,X2,Fext,Options.Gamma,Options.Kappa,Options.Delta);

    if ~resolution_check(X2)
        X2 = reparameterize(X2,Options.nPoints);
        opt_betas = get_betas(X2,X,betas);
        S=SnakeInternalForceMatrix2D_pos(Options.nPoints,Options.Alpha,opt_betas,Options.Gamma);
    end    
    
    % Show current contour
    if(Options.Verbose)
        set(hc2, 'XData', X2(:,2)', 'YData', X2(:,1)');
        drawnow
    end

end

if(nargout>1)
     J=DrawSegmentedArea2D(X2,size(I));
end

end