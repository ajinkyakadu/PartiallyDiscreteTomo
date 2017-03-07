im = imreadgs('cylinders.png');
im = im / max(im(:));

[m,n] = size(im);

nAngles = 15;

proj_geom = astra_create_proj_geom('parallel', 1, max(m,n), linspace2(0,pi,nAngles));
vol_geom  = astra_create_vol_geom(m,n);

W = astra.opTomo('cuda', proj_geom, vol_geom);

p = W*im(:);

sinogram = reshape(p, W.proj_size);
sinogram = astra_add_noise_to_sino_fixed_scaling(sinogram, 5e2);

greyValues = unique(im);

initial_arm_it = 40;
arm_it  = 50;
dart_it = 40;
% This one works the best
mode = 3;
% regularization parameter (is very important!)
lambda = 1;

x_seg = astra.sdart(W, sinogram(:), W.vol_size, greyValues, 'cgls',  ...
    initial_arm_it, arm_it, dart_it, 3, lambda, [], im);

figure
x_dart = astra.dart(sinogram(:), proj_geom, vol_geom, greyValues, initial_arm_it, ...
    arm_it, dart_it, 'SIRT_CUDA', 0.99, [], im);