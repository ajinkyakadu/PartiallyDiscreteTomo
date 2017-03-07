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

% plot sinogram
% figure(10);
% imshow(sinogram,[]);title('clean sinogram')
% figure(11);
% imshow(sinogramN,[]); title('sinogram with noise');
imV = im(:);
imshape = zeros(size(imV));
imshape(imV == 1) = 1;

%% Reconstruction - CGLS

recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
cfg = astra_struct('CGLS_CUDA');
cfg.ProjectionDataId = sinogram_id;
cfg.ReconstructionDataId = recon_id;
cgls_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', cgls_id, 1000);
V = astra_mex_data2d('get', recon_id);
imshow(V, []);

Vvec = V(:);
IV = zeros(size(Vvec));
IV(Vvec >= 1) = 1;
imshow(reshape(IV,256,256),[])

