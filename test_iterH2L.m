% -----------------------------------------------------------------------
% This file is part of the ASTRA Toolbox
% 
% Copyright: 2010-2016, iMinds-Vision Lab, University of Antwerp
%            2014-2016, CWI, Amsterdam
% License: Open Source under GPLv3
% Contact: astra@uantwerpen.be
% Website: http://www.astra-toolbox.com/
% -----------------------------------------------------------------------

% This sample illustrates the use of opTomo.
%
% opTomo is a wrapper around the FP and BP operations of the ASTRA Toolbox,
% to allow you to use them as you would a matrix.
%
% This class requires the Spot Linear-Operator Toolbox to be installed.
% You can download this at http://www.cs.ubc.ca/labs/scl/spot/

startup

% set path to save figures and data
path = strcat(pwd,'/results/limited/');

% set global variables
global fIter
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

% discretize model on fine grid to generate data
n = 256*2;
[im0,~] = createPhantom(0:1/(n-1):1,0:1/(n-1):1,modelOpt); % object of size 256 x 256
x0 = im0(:);

% plot model and background
% fig1 = figure(1); imagesc(im); colorbar; 
% saveas(fig1,'model.fig');
% 
% fig12 = figure(12); imagesc(bgIm); colorbar
% saveas(fig12,'backgroundmodel.fig');

%% Setting up the geometry
% projection geometry
proj_geom = astra_create_proj_geom('parallel', 1, 256, linspace2(0,2*pi/3,5));

% object dimensions
vol_geom  = astra_create_vol_geom(256,256);

% to get rid of inverse crime
% proj_geom0 = astra_create_proj_geom('parallel', 2, 256, linspace2(0,2*pi/3,5));
% vol_geom0  = astra_create_vol_geom(512,512);

%% Generate projection data
% Create the Spot operator for ASTRA using the GPU.
W = opTomo('cuda', proj_geom, vol_geom);

% W0 = opTomo('cuda', proj_geom0, vol_geom0);
p = W*x;

% adding noise to data
pN = addwgn(p,10,0);
% pN = imnoise(p,'poisson');
% pN = imnoise(p,'salt & pepper',0.1);

% reshape the vector into a sinogram
sinogram = reshape(p, W.proj_size);  
sinogramN= reshape(pN, W.proj_size); % look at how noise has been added

% plot sinogram
% figure(10);
% imshow(sinogram,[]);title('clean sinogram')
% figure(11);
% imshow(sinogramN,[]); title('sinogram with noise');

imV = im(:);
imshape = zeros(size(imV));
imshape(imV == 1) = 1;

%% Reconstruction - LSQR
% We use a least squares solver lsqr from Matlab to solve the 
% equation W*x = p.
% Max number of iterations is 100, convergence tolerance of 1e-6.
[y] = lsqr(W, pN, 1e-6, 5000);

rec0shape = zeros(size(imV));
rec0shape(y >= 1) = 1;

% the output is a vector, so we reshape it into an image
reconstruction = reshape(y, W.vol_size);
residual = norm(reconstruction(:) - im(:));
diff0    = rec0shape(:) - imshape(:);
res0Shape= sum(abs(diff0));

fig2 = figure(2);
subplot(2,2,1);
imagesc(reconstruction,[0 1]); axis equal tight; axis off;% imshow(reconstruction, []);
title('Reconstruction');

subplot(2,2,2);
imagesc(im,[0 1]); axis equal tight; axis off;% imshow(im, []);
title('Ground truth');

% The transpose of the operator corresponds to the backprojection.
backProjection = W'*pN;
subplot(2,2,3);
imagesc(reshape(backProjection, W.vol_size)); axis equal tight; axis off;% imshow(reshape(backProjection, W.vol_size), []);
title('Backprojection');

subplot(2,2,4);
imagesc(reshape(diff0,W.vol_size),[-1 1]); axis equal tight; axis off;
title('final difference');

fprintf('\n old residual is : %.2g and old shape residual is %.2g \n',residual,res0Shape);

pause(0.1);

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

fig.show = 1;
fig.save = 1;
fig.path = path;

ipStr.proj_geom = proj_geom;
ipStr.vol_geom  = vol_geom;
ipStr.newMethod = newMethod;
ipStr.fig       = fig;
    
[y1,Op] = jointRec(W,pN,lambda,kappa,maxIter,maxInnerIter,ipStr);

% the output is a vector, so we reshape it into an image
reconstruction1 = reshape(y1, W.vol_size);


rec1shape = zeros(size(imV));
rec1shape(y1 == 1) = 1;

residual1 = norm(reconstruction1(:) - im(:));
difference = rec1shape(:) - imshape(:);
resShape  = sum(abs(difference));

fprintf('\n new residual and resShape is : %0.2g  and  %0.2g \n',residual1,resShape);



fig3 = figure(3);
subplot(2,2,1);
imagesc(reconstruction1,[0 1]); axis equal tight; axis off;% imshow(reconstruction, []);
title('Reconstruction');

subplot(2,2,2);
imagesc(im,[0 1]); axis equal tight; axis off;% imshow(im, []);
title('Ground truth');

% The transpose of the operator corresponds to the backprojection.
backProjection = W'*pN;
subplot(2,2,3);
imagesc(reshape(backProjection, W.vol_size)); axis equal tight; axis off;% imshow(reshape(backProjection, W.vol_size), []);
title('Backprojection');

subplot(2,2,4);
imagesc(reshape(difference,W.vol_size),[-1 1]); axis equal tight; axis off;
title('final difference');


%% saving

if fig.save
    savefig(fig2,strcat(path,'normrec.fig'),'compact');
    saveas(fig2,strcat(path,'normrec'),'epsc');
    savefig(fig3,strcat(path,'jointrec.fig'),'compact');
    saveas(fig3,strcat(path,'jointrec'),'epsc');
end
save(strcat(path,'data.mat'),'y','y1','Op','im','residual','resShape','residual1','lambda');

