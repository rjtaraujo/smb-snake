%%
% This is a script aiming to show the flexibility of the Sparse
% Multi-bending snake. Some examples are provided as a starting point. For
% more information on the method, please read our paper: Ricardo J. Araújo,
% Kelwin Fernandes, Jaime S. Cardoso, "Sparse multi-bending snakes", IEEE
% Transactions on Image Processing. If you have any doubts, please feel free
% to send an email to ricardo.j.araujo@inesctec.pt

%%

close all, clear, clc
addpath('images')
addpath('smb_functions')

% make sure the name below matches the name of the implementation you have 
% from Dirk-Jan Kroon
addpath(genpath('BasicSnake_version2f'))

% select the example (4 are provided)
example = 4;


% these are parameters related to the traditional snake
% General parameters
Options.Verbose = 1;
Options.Iterations = 0; %this will be defined according to each example
Options.nPoints = 0;    %this will be defined according to each example
        
% Edge energy related terms
Options.Sigma1 = 1;
Options.Wedge = 2;
Options.Wline = 0;
Options.Wterm = 0;
Options.Sigma2 = 2;

% GVF related terms
Options.Mu = 0.1;
Options.GIterations = 2000;
Options.Sigma3 = 1;

% Internal energy related terms
Options.Alpha = 0;
Options.Beta = 0;
Options.Delta = 0;      
        

switch example
    case 1
        
        Options.nPoints = 180;
        Options.Iterations = 400;
        Options.betaMean = 5;
        Options.maxValue = 12;
        Options.lambda = 5;
        I = im2double(imread('synth_obj1.png'));

        
        Xi = getInitialContour([100, 90],70,60,Options.nPoints);
        Xf = SMB_Snake2D(I, Xi, Options);
        
    case 2
        
        I = im2double(imread('synth_obj2.png'));
        
        Options.Iterations = 250;
        Options.nPoints = 240;    
        Options.betaMean = 5;
        Options.maxValue = 10;
        Options.lambda = 5;        
        
        Xi = getInitialContour([110, 120],80,90,Options.nPoints);
        Xf = SMB_Snake2D(I, Xi, Options);
        
    case 3

        I = im2double(imread('lung_img.png'));
        
        L = lungPreSegmentation(I);

        Options.Iterations = 250;
        Options.nPoints = 120;    
        Options.betaMean = 10;
        Options.maxValue = 15;
        Options.lambda = 2;        

        
        Xi = getInitialContour([270, 340],120,75,Options.nPoints);
        Xf = SMB_Snake2D(L, Xi, Options);

        
    case 4
        % make sure the name below matches the name of the implementation you have 
        % from Dirk-Jan Kroon
        addpath(genpath('region_growing'))
        
        Options.Iterations = 350;
        Options.nPoints = 220;    
        Options.betaMean = 15;
        Options.maxValue = 40;
        Options.lambda = 20;
        
        I = im2double(imread('hand_depth_img.png'));
        
        H = handRegionGrowing(I);
        
        Xi = getInitialContour([150, 150],70,50,Options.nPoints);        
        Xf = SMB_Snake2D(H, Xi, Options);
        
end
