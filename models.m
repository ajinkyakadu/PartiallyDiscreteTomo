
% model 1
modelOpt.xwidth    = 0.6;
modelOpt.zwidth    = 0.4;
modelOpt.nrand     = 50;
modelOpt.randi     = 6;
modelOpt.bg.smooth = 10;
modelOpt.bg.bmax   = 0.5;
n = 256*2;
[im1,bgIm1] = createPhantom(0:1/(n-1):1,0:1/(n-1):1,modelOpt); % object of size 256 x 256

figure(1);imagesc(im1,[0 1]); axis equal tight; title('Model')
figure(2);imagesc(bgIm1,[0 1]); axis equal tight;title('Background');

%% model 2
modelOpt.xwidth    = 0.6;
modelOpt.zwidth    = 0.4;
modelOpt.nrand     = 35;
modelOpt.randi     = 90;
modelOpt.bg.smooth = 5;
modelOpt.bg.bmax   = 0.5;
n = 256*2;
[im2,bgIm2] = createPhantom(0:1/(n-1):1,0:1/(n-1):1,modelOpt);

figure(3);imagesc(im2,[0 1]); axis equal tight; title('Model')
figure(4);imagesc(bgIm2,[0 1]); axis equal tight;title('Background');

%% model 3
modelOpt.xwidth    = 0.5;
modelOpt.zwidth    = 0.3;
modelOpt.nrand     = 43;
modelOpt.randi     = 56;
modelOpt.bg.smooth = 20;
modelOpt.bg.bmax   = 0.5;
n = 256*2;
[im3,bgIm3] = createPhantom(0:1/(n-1):1,0:1/(n-1):1,modelOpt);

figure(6);imagesc(im3,[0 1]); axis equal tight; title('Model')
% figure(6);imagesc(bgIm3,[0 1]); axis equal tight;title('Background');

%% model 4
[x,y] = ginput(4);


% model 5