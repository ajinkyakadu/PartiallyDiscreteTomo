%% 

startup

% set path to save figures and data
path = strcat(pwd,'/results/phantom3/nm/');

% set global variables
global fIter
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
modelOpt.type      = 3;
n = 128;
[im,bgIm] = createPhantom(0:1/(n-1):1,0:1/(n-1):1,modelOpt); % object of size 256 x 256

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

% reshape the vector into a sinogram
sinogram = reshape(pN, W.proj_size);  

imV = im(:);
imshape = zeros(size(imV));
imshape(imV == 1) = 1;

%% Reconstruction - LSQR
[x_ls] = lsqr(W, pN, 1e-6, 5000);

LS.shape = zeros(size(x_ls));
LS.shape(x_ls >= 1) = 1;

LS.modRes   = norm(x_ls - im(:));
LS.diff     = LS.shape - imshape;
LS.shapeRes = sum(abs(LS.diff));
LS.dataRes  = norm(W*x_ls - pN);

fprintf('\n LSQR Method: ModelResidual = %0.2g and ShapeResidual =  %0.2g DataResidual = %0.2g \n',LS.modRes,LS.shapeRes,LS.dataRes);


fig1 = plotRec(x_ls,W,pN,im,LS.diff);


%% Reconstruction with new method
% We use a joint reconstruction method to solve the 
% equation W*x = p.
% jointRec(W,p,lambda,kappa,maxIter,iPstr);


lambda      = 2.976e4;
kappa       = 0.05;
maxIter     = 50;
maxInnerIter= 100;

newMethod.run       = 0;
newMethod.maxLoop   = 3;
newMethod.changeL   = 0;
newMethod.changeHW  = 1;
newMethod.maxIter   = 20;
newMethod.recFactL  = 0.1;
newMethod.recFactHW = 0.9;

ipStr.proj_geom = proj_geom;
ipStr.vol_geom  = vol_geom;
ipStr.newMethod = newMethod;
ipStr.fig       = fig;
    
[x_nm,Op] = jointRec(W,pN,lambda,kappa,maxIter,maxInnerIter,ipStr);

NM.shape = zeros(size(x_nm));
NM.shape(x_nm >= 1) = 1;

NM.modRes   = norm(x_nm - im(:));
NM.diff     = NM.shape - imshape;
NM.shapeRes = sum(abs(NM.diff));
NM.dataRes  = norm(W*x_nm - pN);

fprintf('\n New Method: ModelResidual = %0.2g and ShapeResidual =  %0.2g DataResidual = %0.2g \n',NM.modRes,NM.shapeRes,NM.dataRes);

fig2 = plotRec(x_nm,W,pN,im,NM.diff);



%% saving

if fig.save
    savefig(fig1,strcat(path,'LSQR.fig'),'compact');
    saveas(fig1,strcat(path,'LSQR'),'epsc');
    savefig(fig2,strcat(path,'NewMethod.fig'),'compact');
    saveas(fig2,strcat(path,'NewMethod'),'epsc');
end
save(strcat(path,'data.mat'),'x_ls','x_nm','LS','NM');

