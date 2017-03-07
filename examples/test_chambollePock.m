%% setup
n = 512;
nAngles = 20;

im = phantom(n);

proj_geom = astra_create_proj_geom('parallel', 1, n, linspace2(0,pi,nAngles));
vol_geom  = astra_create_vol_geom(size(im));

W = astra.opTomo('cuda', proj_geom, vol_geom);

p = W*im(:);

TV = astra.tv.opTV(n);

%% reconstruct
scale = 1/30;
x_tv = astra.tv.chambollePock(scale*W, TV, p, 200, 1e-1, true, [], 2);

%% least squares
x_ls = lsqr(W,p);

%% visualize
subplot(1,2,1);
show(x_ls);

subplot(1,2,2);
show(x_tv);
