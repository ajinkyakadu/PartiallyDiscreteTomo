

clearvars; close all;

global fIter 
fIter   = 1;
mtype   = 1;
noise   = 1;

path    = strcat(pwd,'/results/');

%% load a phantom image
if (mtype==1)
    n   = 256;  mtype1  = 1;    rseed   = 1;    bgmax   = 0.5;
    lambdaPLS   = 2.976e7;    lambdaTV  = 0.1129;
    greyValues      = [linspace(0,0.5,20) 1]';
elseif (mtype==2)
    n   = 128;  mtype1  = 2;    rseed   = 5;    bgmax   = 0.5;
    lambdaPLS   = 1.129e9;    lambdaTV  = 0.1129;
    greyValues      = [linspace(0,0.5,20) 1]';
elseif (mtype==3)
    n   = 256;  mtype1  = 1;    rseed   = 10;   bgmax   = 0.8;
    lambdaPLS   = 2.976e7;    lambdaTV  = 0.1129;
    greyValues      = [linspace(0,0.8,30) 1]';
elseif (mtype==4)
    n   = 128;  mtype1  = 2;    rseed   = 20;   bgmax   = 0.8;
    lambdaPLS   = 1.274e8;    lambdaTV  = 0.04833;
    greyValues      = [linspace(0,0.8,30) 1]';
elseif (mtype==5)
    n   = 128;  mtype1  = 3;    rseed   = 30;   bgmax   = 0;
elseif (mtype==6)
    n   = 256;  mtype1  = 1;    rseed   = 30;   bgmax   = 0;
end


modelOpt.xwidth     = 0.6;
modelOpt.zwidth     = 0.4;
modelOpt.nrand      = 50;
modelOpt.randi      = 6;
modelOpt.bg.smooth  = 10;
modelOpt.bg.bmax    = bgmax;
modelOpt.type       = mtype1;
modelOpt.rseed      = rseed;
modelOpt.gV         = greyValues(1:end-1);

[im,bgIm]   = createPhantom(0:1/(n-1):1,0:1/(n-1):1,modelOpt); % object of size 256 x 256
x           = im(:);

fig1 = figure(1); imagesc(im,[0 1]); axis equal tight; axis off; colormap gray
saveas(fig1,strcat(path,'model',num2str(mtype)),'epsc');
saveas(fig1,strcat(path,'model',num2str(mtype)),'fig');
pause(0.001);

imV = im(:);
imshape = zeros(size(imV));
imshape(imV == 1) = 1;

%% Setting up the geometry
% projection geometry
proj_geom = astra_create_proj_geom('parallel', 1, n, linspace2(0,pi,180));

% object dimensions
vol_geom  = astra_create_vol_geom(n,n);


%% Generate projection data
% Create the Spot operator for ASTRA using the GPU.
W   = opTomo('cuda', proj_geom, vol_geom);

W0  = opTomo('line', proj_geom, vol_geom);
p   = W0*x;

% adding noise to data
if noise
    pN = addwgn(p,3,0);
else
    pN = p;
end

% reshape the vector into a sinogram
sinogram = reshape(p, W.proj_size);  
sinogramN= reshape(pN, W.proj_size); % look at how noise has been added

%% Reconstruction - LSQR
% We use a least squares solver lsqr from Matlab to solve the 
% equation W*x = p.
% Max number of iterations is 100, convergence tolerance of 1e-6.
% [x_ls]  = lsqr(W, pN, 1e-6, 5000);
% rec_ls  = reshape(x_ls, W.vol_size);
% res_ls  = norm(rec_ls(:) - im(:));
% 
% fig2 = figure(2);
% imagesc(rec_ls,[0 1]); axis equal tight; axis off;% imshow(reconstruction, []);
% saveas(fig2,strcat(path,'m',num2str(mtype),'_lsqr_n',num2str(noise)),'epsc');
% saveas(fig2,strcat(path,'m',num2str(mtype),'_lsqr_n',num2str(noise)),'fig')
% 
% LS.rec              = rec_ls;
% LS.shape            = zeros(size(x_ls));
% LS.shape(x_ls >= 1) = 1;
% LS.modRes           = norm(x_ls - im(:));
% LS.diff             = LS.shape - imshape;
% LS.shapeRes         = sum(abs(LS.diff));
% LS.dataRes          = norm(W*x_ls - pN);
% 
% fprintf('\n LSQR Method: ModelResidual = %0.2d and ShapeResidual =  %0.2d DataResidual = %0.2d \n',LS.modRes,LS.shapeRes,LS.dataRes);

%% Reconstruction with parametric level-set method
% We use a joint reconstruction method to solve the 
% equation W*x = p.
% jointRec(W,p,lambda,kappa,maxIter,iPstr)
kappa               = 0.05/8;
maxIter             = 50;
maxInnerIter        = 200;
fig.show            = 0;
fig.save            = 0;
fig.path            = path;
ipStr.fig           = fig;


% RBF Kernel
Koptions.tau        = 5;        % how coarse the RBF grid should be wrt computational grid
Koptions.eta        = 4;        % parameter to control the spread of RBF
Koptions.nouter     = 2;        % RBF layers outside compuational domain 
Koptions.rtype      = 'compact';% RBF type
Koptions.ltype      = 'L2';     % distance norms for RBF
xt                  = 0:(1/(n-1)):1;    % x-dn vector 
zt                  = 0:(1/(n-1)):1;    % z-dn vector
[zz,xx]             = ndgrid(zt,xt);
[A,nr,Xc,Zc]        = generateKernel(xt,zt,Koptions);   
A                   = opMatrix(A);

% initialize alpha, level-set parameter
a0                  = -1*ones(nr);
a0(floor(nr(1)/2)-3:floor(nr(1)/2),floor(nr(2)/2):floor(nr(2)/2)+3) = 1;    % initial level-set
a0                  = a0(:);

% regularizer
L                   = opTV(n);

% fig11 = figure(11);surfl(xx,zz,reshape(A*a0,[n n]));shading interp; hold on;
% hfig1 = imagesc(xt,zt,reshape(A*a0,[n n])); hfig1.AlphaData = 0.25; axis([0 1 0 1])
% axis ij; axis off;
% text(0,-0.5,-2,'x','HorizontalAlignment','left','FontSize',12);
% text(-0.5,0,-2,'y','HorizontalAlignment','left','FontSize',12); hold off;
% saveas(fig11,strcat(path,'initLSF'),'png');
% saveas(fig11,strcat(path,'initLSF'),'fig');

% Xcp = Xc(a0>0);
% Zcp = Zc(a0>0);
% Xcn = Xc(a0<0);
% Zcn = Zc(a0<0);
% fig12 = figure(12);
% plot(Xcp,Zcp,'r+','MarkerSize',4);hold on;
% plot(Xcn,Zcn,'bo','MarkerSize',2); axis ij;
% contour(xx,zz,reshape(A*a0,[n n]),[0 0],'LineWidth',2); hold off;
% xlabel('x','FontSize',12); ylabel('y','FontSize',12);
% saveas(fig12,strcat(path,'initLS'),'epsc');
% saveas(fig12,strcat(path,'initLS'),'fig');

[x_pls,Op]  = jointRec(W,pN,a0,A,L,lambdaPLS,kappa,maxIter,maxInnerIter,ipStr);
rec_pls     = reshape(x_pls, W.vol_size);
af          = Op.xf;

fig3 = figure(3);
imagesc(rec_pls,[0 1]); axis equal tight; axis off; colormap gray
saveas(fig3,strcat(path,'m',num2str(mtype),'_pls_n',num2str(noise)),'epsc');
saveas(fig3,strcat(path,'m',num2str(mtype),'_pls_n',num2str(noise)),'fig')

% fig13 = figure(13);surfl(xx,zz,reshape(A*af,[n n]));shading interp; hold on;
% hfig2 = imagesc(xt,zt,reshape(A*af,[n n])); hfig2.AlphaData = 0.25; axis([0 1 0 1])
% axis ij; axis off;
% text(0,-0.5,-2,'x','HorizontalAlignment','left','FontSize',12);
% text(-0.5,0,-2,'y','HorizontalAlignment','left','FontSize',12); hold off;
% saveas(fig13,strcat(path,'finalLSF'),'png');
% saveas(fig13,strcat(path,'finalLSF'),'fig');

% Xcp1 = Xc(af>0);
% Zcp1 = Zc(af>0);
% Xcn1 = Xc(af<0);
% Zcn1 = Zc(af<0);
% fig14 = figure(14);
% plot(Xcp1,Zcp1,'r+','MarkerSize',4);hold on;
% plot(Xcn1,Zcn1,'bo','MarkerSize',2);
% contour(xx,zz,reshape(A*af,[n n]),[0 0],'LineWidth',2);hold off;
% xlabel('x','FontSize',12); ylabel('y','FontSize',12);
% saveas(fig14,strcat(path,'finalLS'),'epsc');
% saveas(fig14,strcat(path,'finalLS'),'fig');

PLS.rec                 = rec_pls;
PLS.shape               = zeros(size(x_pls));
PLS.shape(x_pls >= 1)   = 1;
PLS.modRes              = norm(x_pls - im(:));
PLS.diff                = PLS.shape - imshape;
PLS.shapeRes            = sum(abs(PLS.diff));
PLS.dataRes             = norm(W*x_pls - pN);
PLS.jacId               = nnz(PLS.shape & imshape)/nnz(PLS.shape | imshape);

fprintf('\n PLS Method: ModelResidual = %0.2d and ShapeResidual =  %0.2d DataResidual = %0.2d Jaccard Index = %0.2d \n',PLS.modRes,PLS.shapeRes,PLS.dataRes,PLS.jacId);


%% Total Variation
TVOp    = opTV(n);
scale   = sqrt(eigs(TVOp'*TVOp,1)/eigs(W'*W,1));


x_tv    = chambollePock(scale*W, TVOp, pN, 200, lambdaTV, true, [], 0);   % lambda = 1.9 best MR, = 2.683 best SR
x_tv    = scale*x_tv;
rec_tv  = reshape(x_tv, W.vol_size);

for i=0.6:0.01:1
    tvshape = zeros(size(x_tv));
    tvshape(x_tv>=i) = 1;
    p = int8(100*(i-0.6)+1);
    jacId(p) = nnz(tvshape & imshape)/nnz(tvshape | imshape);
end

[~,j_mid] = max(jacId);
thr = 0.01*(j_mid-1)+0.6;

fig4 = figure(4);
imagesc(rec_tv,[0 1]); axis equal tight; axis off; colormap gray
hold on; contour(rec_tv,[thr thr],':r','LineWidth',2); hold off;
saveas(fig4,strcat(path,'m',num2str(mtype),'_tv_n',num2str(noise)),'epsc');
saveas(fig4,strcat(path,'m',num2str(mtype),'_tv_n',num2str(noise)),'fig')

TV.rec              = rec_tv;
TV.shape            = zeros(size(x_tv));
TV.shape(x_tv >=thr)= 1;
TV.modRes           = norm(x_tv - im(:));
TV.diff             = TV.shape - imshape;
TV.shapeRes         = sum(abs(TV.diff));
TV.dataRes          = norm(W*x_tv - pN);
TV.jacId            = nnz(TV.shape & imshape)/nnz(TV.shape | imshape);

fprintf('\n Total Variation Method: ModelResidual = %0.2d ShapeResidual =  %0.2d DataResidual = %0.2d Jaccard Index = %0.2d \n',TV.modRes,TV.shapeRes,TV.dataRes,TV.jacId);

%% DART


initial_arm_it  = 40;
arm_it          = 5;
dart_it         = 40;

x_dart      = astra.dart(sinogramN(:), proj_geom, vol_geom, greyValues, initial_arm_it, ...
            arm_it, dart_it, 'SIRT_CUDA', 0.99, [], im);
x_dart      = x_dart(:);
rec_dart    = reshape(x_dart, W.vol_size);

fig5 = figure(5);
imagesc(rec_dart,[0 1]); axis equal tight; axis off; colormap gray
saveas(fig5,strcat(path,'m',num2str(mtype),'_dart_n',num2str(noise)),'epsc');
saveas(fig5,strcat(path,'m',num2str(mtype),'_dart_n',num2str(noise)),'fig')

DART.rec                = rec_dart;
DART.shape              = zeros(size(x_dart));
DART.shape(x_dart >= 1) = 1;
DART.modRes             = norm(x_dart - im(:));
DART.diff               = DART.shape - imshape;
DART.shapeRes           = sum(abs(DART.diff));
DART.dataRes            = norm(W*x_dart - pN);
DART.jacId              = nnz(DART.shape & imshape)/nnz(DART.shape | imshape);

fprintf('\n DART Method: ModelResidual = %0.2d and ShapeResidual =  %0.2d DataResidual = %0.2d Jaccard Index = %0.2d \n',DART.modRes,DART.shapeRes,DART.dataRes,DART.jacId);

%% P-DART

mask_id = astra_mex_data2d('create','-vol',vol_geom,1);

[sinogram_id, sinogram] = astra_create_sino_cuda(im, proj_geom, vol_geom);
% proj_id = astra_mex_data2d('create','-vol',proj_geom,p);
vol_id = astra_mex_data2d('create','-vol',vol_geom,0);

if mtype<3; tau=0.5; else; tau = 0.8; end;
rho = 1;

cfg = astra_struct('SIRT_CUDA');
cfg.ProjectionDataId = sinogram_id;
cfg.ReconstructionDataId = vol_id;
cfg.option.ReconstructionMaskId = mask_id;
alg_id = astra_mex_algorithm('create',cfg);

for i=1:150
    astra_mex_algorithm('iterate',alg_id,1);
    x_pdart = astra_mex_data2d('get',vol_id);
    
    q = x_pdart < tau;
    astra_mex_data2d('set',mask_id,q);
    
    segmentation = double(~q)*rho;
    x_pdart(~q) = rho;
    astra_mex_data2d('set',vol_id,x_pdart);
    
    [tmp_id,Ws] = astra_create_sino_cuda(segmentation,proj_geom,vol_geom);
    astra_mex_data2d('set',sinogram_id,sinogram-Ws);
end

x_pdart     = x_pdart(:);
rec_pdart   = reshape(x_pdart, W.vol_size);

fig6 = figure(6);
imagesc(rec_pdart,[0 1]); axis equal tight; axis off; colormap gray
saveas(fig6,strcat(path,'m',num2str(mtype),'_pdart_n',num2str(noise)),'epsc');
saveas(fig6,strcat(path,'m',num2str(mtype),'_pdart_n',num2str(noise)),'fig')

PDART.rec                   = rec_pdart;
PDART.shape                 = zeros(size(x_pdart));
PDART.shape(x_dart >= 1)    = 1;
PDART.modRes                = norm(x_pdart - im(:));
PDART.diff                  = PDART.shape - imshape;
PDART.shapeRes              = sum(abs(PDART.diff));
PDART.dataRes               = norm(W*x_pdart - pN);
PDART.jacId                 = nnz(PDART.shape & imshape)/nnz(PDART.shape | imshape);


fprintf('\n P-DART Method: ModelResidual = %0.2d and ShapeResidual =  %0.2d DataResidual = %0.2d Jaccard Index = %0.2d \n',PDART.modRes,PDART.shapeRes,PDART.dataRes,PDART.jacId);

%% saving

save(strcat(path,'results',num2str(noise),'.mat'),'PLS','TV','DART','PDART');


