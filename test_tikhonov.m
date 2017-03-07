
startup

global fIter 
path = strcat(pwd,'/results/tikhonov/');
fIter = 1;

fig.show = 1;
fig.save = 1;
fig.path = path;


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

x = im(:);


%% Setting up the geometry
% projection geometry
proj_geom = astra_create_proj_geom('parallel', 1, 256, linspace2(0,2*pi/3,5));

% object dimensions
vol_geom  = astra_create_vol_geom(256,256);

%% Generate projection data
% Create the Spot operator for ASTRA using the GPU.
W = opTomo('cuda', proj_geom, vol_geom);

p = W*x;

% adding noise to data
pN = addwgn(p,10,0);

% reshape the vector into a sinogram
sinogram = reshape(pN, W.proj_size);  

imV = im(:);
imshape = zeros(size(imV));
imshape(imV == 1) = 1;


%% Reconstruction - LSQR
[x_ls] = lsqr(W, pN, 1e-6, 500);

LS.shape = zeros(size(x_ls));
LS.shape(x_ls >= 1) = 1;

LS.modRes   = norm(x_ls - im(:));
LS.diff     = LS.shape - imshape;
LS.shapeRes = sum(abs(LS.diff));
LS.dataRes  = norm(W*x_ls - pN);

fprintf('\n LSQR Method: ModelResidual = %0.2g and ShapeResidual =  %0.2g DataResidual = %0.2g \n',LS.modRes,LS.shapeRes,LS.dataRes);

fig1 = plotRec(x_ls,W,pN,im,LS.diff);

%% Reconstruction with Tikhonov Regularization


regMode = 1;    % type of regularizer/penalty, 1: Tikhonov for smoothness, 2: penalty for x - x0
L = smoothReg([256 256],regMode);


fig.show = 0;
fig.save = 1;

lambda = 10; % logspace(-10,10,42);
H1 = speye(256*256);
L1 = sqrt(lambda)*L;

Wu = opTomoU('cuda', proj_geom, vol_geom,H1,L1);
[x_tr] = lsqr(Wu, [pN;zeros(256*256,1)], 1e-6, 500);

TR.shape = zeros(size(x_tr));
TR.shape(x_tr >= 1) = 1;

TR.modRes   = norm(x_tr - im(:));
TR.diff     = TR.shape - imshape;
TR.shapeRes = sum(abs(TR.diff));
TR.dataRes  = norm(W*x_tr - pN);

fprintf('\n Tikhonov Reg Method: ModelResidual = %0.2g and ShapeResidual =  %0.2g DataResidual = %0.2g \n',TR.modRes,TR.shapeRes,TR.dataRes);

fig2 = plotRec(x_tr,W,pN,im,LS.diff);

%% saving

if fig.save
    savefig(fig1,strcat(path,'LSQR.fig'),'compact');
    saveas(fig1,strcat(path,'LSQR'),'epsc');
    savefig(fig2,strcat(path,'TR.fig'),'compact');
    saveas(fig2,strcat(path,'TR'),'epsc');
end
save(strcat(path,'data.mat'),'x_ls','x_tr','LS','TR');

