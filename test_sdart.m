startup

global fIter 
path = strcat(pwd,'/results/sdart/');
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
% and flatten it to a vector
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


imV = im(:);
imshape = zeros(size(imV));
imshape(imV == 1) = 1;

sinogram = reshape(pN, W.proj_size);
% sinogram = astra_add_noise_to_sino_fixed_scaling(sinogram, 5e2);

greyValues = unique(im);

%% reconstruction - SDART

initial_arm_it = 40;
arm_it  = 50;
dart_it = 40;
% This one works the best
mode = 3;
% regularization parameter (is very important!)
lambda = 1;

x_sdart = astra.sdart(W, sinogram(:), W.vol_size, greyValues, 'lsqr',  ...
    initial_arm_it, arm_it, dart_it, 3, lambda, [], im);

SD.shape = zeros(size(x_sdart));
SD.shape(x_sdart >= 1) = 1;

SD.modRes   = norm(x_sdart - im(:));
SD.diff     = SD.shape - imshape;
SD.shapeRes = sum(abs(SD.diff));
SD.dataRes  = norm(W*x_sdart - pN);

fprintf('\n SDART Method: ModelResidual = %0.2g and ShapeResidual =  %0.2g DataResidual = %0.2g \n',SD.modRes,SD.shapeRes,SD.dataRes);


fig1 = plotRec(x_sdart,W,pN,im,SD.diff);

%% reconstruction - DART
figure
x_dart = astra.dart(sinogram(:), proj_geom, vol_geom, greyValues, initial_arm_it, ...
   arm_it, dart_it, 'SIRT_CUDA', 0.99, [], im);
x_dart = x_dart(:);

DT.shape = zeros(size(x_dart));
DT.shape(x_dart >= 1) = 1;

DT.modRes   = norm(x_dart - im(:));
DT.diff     = DT.shape - imshape;
DT.shapeRes = sum(abs(DT.diff));
DT.dataRes  = norm(W*x_dart - pN);

fprintf('\n DART Method: ModelResidual = %0.2g and ShapeResidual =  %0.2g DataResidual = %0.2g \n',DT.modRes,DT.shapeRes,DT.dataRes);


fig2 = plotRec(x_dart,W,pN,im,DT.diff);

%% saving

if fig.save
    savefig(fig1,strcat(path,'SDART'),'compact');
    saveas(fig1,strcat(path,'SDART'),'epsc');
    savefig(fig2,strcat(path,'DART'),'compact');
    saveas(fig2,strcat(path,'DART'),'epsc');
end
save(strcat(path,'data.mat'),'x_sdart','x_dart','SD','DT');