

startup

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

%% create forward projection
[sinogram_id, sinogram] = astra_create_sino_cuda(im, proj_geom, vol_geom);

%% reconstruct
recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
cfg = astra_struct('SART_CUDA');
cfg.ProjectionDataId = sinogram_id;
cfg.ReconstructionDataId = recon_id;
cfg.option.ProjectionOrder = 'custom';
cfg.option.ProjectionOrderList = [0:5:175 1:5:176 2:5:177 3:5:178 4:5:179];
sart_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', sart_id, 500);
V = astra_mex_data2d('get', recon_id);
imshow(V, []);