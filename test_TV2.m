startup

global fIter 
path = strcat(pwd,'/results/phantom2/TV2/');
fIter = 1;

fig.show = 0;
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
modelOpt.type      = 2;
n = 128;
[im,bgIm] = createPhantom(0:1/(n-1):1,0:1/(n-1):1,modelOpt); % object of size 256 x 256
% and flatten it to a vector
x = im(:);



%% Setting up the geometry
% projection geometry
proj_geom = astra_create_proj_geom('parallel', 1, n, linspace2(0,2*pi/3,5));

% object dimensions
vol_geom  = astra_create_vol_geom(n,n);


%% Generate projection data
% Create the Spot operator for ASTRA using the GPU.
W = opTomo('cuda', proj_geom, vol_geom);

p = W*x;

% adding noise to data
pN = addwgn(p,10,0);


imV = im(:);
imshape = zeros(size(imV));
imshape(imV == 1) = 1;



%% reconstruct - Total Variation

TVOp = opTV(n);

scale = 0.0796;
x_tv = chambollePock(scale*W, TVOp, pN, 200, 1.9, true, [], 1);   % lambda = 1.9 best MR, = 2.683 best SR
x_tv = scale*x_tv;

%% 

TV.shape = zeros(size(x_tv));
TV.shape(x_tv >= 1) = 1;

TV.modRes   = norm(x_tv - im(:));
TV.diff     = TV.shape - imshape;
TV.shapeRes = sum(abs(TV.diff));
TV.dataRes  = norm(W*x_tv - pN);

fprintf('\n Total Variation Method: ModelResidual = %0.2g and ShapeResidual =  %0.2g DataResidual = %0.2g \n',TV.modRes,TV.shapeRes,TV.dataRes);


fig2 = plotRec(x_tv,W,pN,im,TV.diff);

%% least squares
x_ls = lsqr(W,pN,1e-6,500);

LS.shape = zeros(size(x_ls));
LS.shape(x_ls >= 1) = 1;

LS.modRes   = norm(x_ls - im(:));
LS.diff     = LS.shape - imshape;
LS.shapeRes = sum(abs(LS.diff));
LS.dataRes  = norm(W*x_ls - pN);

fprintf('\n LSQR Method: ModelResidual = %0.2g and ShapeResidual =  %0.2g DataResidual = %0.2g \n',LS.modRes,LS.shapeRes,LS.dataRes);


fig1 = plotRec(x_ls,W,pN,im,LS.diff);

%% saving

if fig.save
    savefig(fig1,strcat(path,'LSQR.fig'),'compact');
    saveas(fig1,strcat(path,'LSQR'),'epsc');
    savefig(fig2,strcat(path,'TV1.fig'),'compact');
    saveas(fig2,strcat(path,'TV1'),'epsc');
end
save(strcat(path,'data1.mat'),'x_ls','x_tv','LS','TV');

