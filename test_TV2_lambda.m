startup

global fIter 
path = strcat(pwd,'/results/tv_lambda/');
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

TVOp = opTV(n);

%% reconstruct
lambda      = logspace(-5,2,50);
TV.modRes    = zeros(size(lambda));
TV.shapeRes  = zeros(size(lambda));
TV.dataRes   = zeros(size(lambda));

for i=1:length(lambda)
    scale = 0.0796;
    x_tv = chambollePock(scale*W, TVOp, pN, 200, lambda(i), true, [], 0);
    x_tv = scale*x_tv;
    
    %%
    
    TV.shape = zeros(size(x_tv));
    TV.shape(x_tv >= 1) = 1;
    
    TV.modRes(i)     = norm(x_tv - im(:));
    TV.diff          = TV.shape - imshape;
    TV.shapeRes(i)   = sum(abs(TV.diff));
    TV.dataRes(i)    = norm(W*x_tv - pN);
    
    fprintf('\n TV Method: lambda %.2g ModelResidual = %.2g and ShapeResidual =  %.2g DataResidual = %.2g \n',lambda(i),TV.modRes(i),TV.shapeRes(i),TV.dataRes(i));
end


%% visualize
fig1 = figure(1);subplot(3,1,1);semilogx(lambda,TV.modRes);
xlabel('\lambda');ylabel('Model Residual');
subplot(3,1,2);semilogx(lambda,TV.shapeRes);
xlabel('\lambda');ylabel('Shape Residual');
subplot(3,1,3);loglog(lambda,TV.dataRes);
xlabel('\lambda');ylabel('Data Residual');
savefig(fig1,strcat(path,'residuals')); saveas(fig1,strcat(path,'residuals'),'epsc');

save('data.mat','lambda','TV');
