%% 

startup

% set path to save figures and data
path = strcat(pwd,'/results/phantom3/NM_lambda/');

% set global variables
global fIter
fIter = 1;

fig.show = 0;
fig.save = 0;
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


%% Reconstruction with new method
% We use a joint reconstruction method to solve the 
% equation W*x = p.
% jointRec(W,p,lambda,kappa,maxIter,iPstr);


lambda      = 1e6;
kappa       = 0.05;
maxIter     = 50;
maxInnerIter= 100;

newMethod.run       = 0;
newMethod.maxLoop   = 10;
newMethod.changeL   = 1;
newMethod.changeHW  = 1;
newMethod.maxIter   = 20;
newMethod.recFactL  = 0.1;
newMethod.recFactHW = 0.9;

ipStr.proj_geom = proj_geom;
ipStr.vol_geom  = vol_geom;
ipStr.newMethod = newMethod;
ipStr.fig       = fig;

lambda      = logspace(3,7,20);
NM.modRes    = zeros(size(lambda));
NM.shapeRes  = zeros(size(lambda));
NM.dataRes   = zeros(size(lambda));

for i=1:length(lambda)
    [x_nm,Op] = jointRec(W,pN,lambda(i),kappa,maxIter,maxInnerIter,ipStr);

    NM.shape = zeros(size(x_nm));
    NM.shape(x_nm >= 1) = 1;

    NM.modRes(i)   = norm(x_nm - im(:));
    NM.diff        = NM.shape - imshape;
    NM.shapeRes(i) = sum(abs(NM.diff));
    NM.dataRes(i)  = norm(W*x_nm - pN);

    fprintf('\n New Method: lambda %0.2g ModelResidual = %0.2g and ShapeResidual =  %0.2g DataResidual = %0.2g \n',lambda(i),NM.modRes(i),NM.shapeRes(i),NM.dataRes(i));
end



%% visualize

fig1 = figure(1);subplot(3,1,1);semilogx(lambda,NM.modRes);
xlabel('\lambda');ylabel('Model Residual');
subplot(3,1,2);semilogx(lambda,NM.shapeRes);
xlabel('\lambda');ylabel('Shape Residual');
subplot(3,1,3);loglog(lambda,NM.dataRes);
xlabel('\lambda');ylabel('Data Residual');
savefig(fig1,strcat(path,'residuals')); saveas(fig1,strcat(path,'residuals'),'epsc');

save('data.mat','lambda','NM');

