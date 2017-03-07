

startup

% addpath(genpath('/home/ajinkya/TomographyJointRC/TVReg/'));
global fIter 
path = strcat(pwd,'/results/tikhonov/l10/');
fIter = 1;


%% load a phantom image

% im = phantom(256);
modelOpt.xwidth    = 0.6;
modelOpt.zwidth    = 0.4;
modelOpt.nrand     = 50;
modelOpt.randi     = 6;
modelOpt.bg.smooth = 10;
modelOpt.bg.bmax   = 0.5;
n = 256;
[im,bgIm] = createPhantom(0:1/(n-1):1,0:1/(n-1):1,modelOpt); % object of size 256 x 256
% and flatten it to a vector
x = im(:);


%% Setting up the geometry
% projection geometry
proj_geom = astra_create_proj_geom('parallel', 1, 256, linspace2(0,2*pi/3,5));

% object dimensions
vol_geom  = astra_create_vol_geom(256,256);

x0 = zeros(size(x));

%% Generate projection data
% Create the Spot operator for ASTRA using the GPU.
W = opTomo('cuda', proj_geom, vol_geom);

p = W*x;

% adding noise to data
pN = addwgn(p,10,0);


%% Set parameters
rnl    = 0.01;       % Noise level
alpha  = 0.1;        % Regularization parameter
r1_max = 13;         % Halfwidth of object cube
N      = 2*r1_max+1; % Full width of object cube
N3     = N^3;        % Total number of variables
dims   = [256 256];    % Dimensions
u_max  = 15;         % Halfwidth of projection planes
U      = 2*u_max+1;  % Full width of projection planes
numProjections = 25; % Number of projections to reconstruct from

%% Parameters for the reconstruction algorithms

tau         = 1e-4*norm(x0,'inf');       % Huber smoothing parameter

% Specify nonnegativity constraints
constraint.type = 2;
constraint.c    = 0*ones(prod(dims),1);
constraint.d    = 1*ones(prod(dims),1);

% Options
opt.epsb_rel = 1e-6;
opt.k_max    = 10000;
opt.qs       = 1;
opt.K        = 2;
opt.verbose  = 1;
opt.beta     = 0.95;

% Options for reference solution
opt_ref.epsb_rel = 1e-8;
opt_ref.k_max    = 20000;
opt_ref.verbose  = 1;

%% Solve: Compute TV minimizer

% Reference solution
[x_ref fxk_ref hxk_ref gxk_ref fxkl_ref info_ref] = ...
    tvreg_upn(W,pN,alpha,tau,dims,constraint,opt_ref);
fs = fxkl_ref(end);     % Final reference objective function value

% Solve using GPBB
tic
[xk_GPBB fxk_GPBB hxk_GPBB gx_kGPBB fxkl_GPBB info_GPBB] = ...
    tvreg_gpbb(A,b,alpha,tau,dims,constraint,opt);
tGPBB = toc

% Solve using UPN
tic
[xk_UPN fxk_UPN hxk_UPN gxk_UPN fxkl_UPN info_UPN] = ...
    tvreg_upn(A,b,alpha,tau,dims,constraint,opt);
tupn = toc

% Solve using UPN0
tic
opt.qs = 0;
[xk_UPNz fxk_UPNz hxk_UPNz gxk_UPNz fxkl_UPNz info_UPNz] = ...
    tvreg_upn(A,b,alpha,tau,dims,constraint,opt);
tupnz = toc