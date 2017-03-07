%% Example 1
% 
%   This example shows reconstruction on the 'optimal' QT grid
%
%

%% phantom
I = imread('../phantoms/spiral_128.png');
n = size(I,1);
I = double(I(:)/max(I(:)));


%% get projections
nr     = n;
ntheta = 32;
theta  = linspace(0,pi,ntheta);

%% paramters
a      = 0;    % noise level
sigma  = 1e-6; % stopping criterion
niter  = 100;
method = 'LSQR';
omega  = .5;

%% get projections
proj_geom = astra_create_proj_geom('parallel', 1.0, n, theta);
vol_geom  = astra_create_vol_geom(n,n);
W = opTomo('strip', proj_geom, vol_geom);
p = W*I(:);
noise = randn(size(p));
p = p + a*noise/norm(noise)*norm(p);

%% fine-scale inversion
If = smoother(W,p,zeros(n^2,1),niter,sigma,omega,method);

%% QT inversion
% get optimal grid
[~, Sk] = qtget(reshape(I,n,n),0, [1 128]);
Vk      = opFunction(n*n,nnz(Sk),@(x,mode)qtProlong(x,Sk,mode));
xk      = smoother(W*Vk,p,Vk'*zeros(n^2,1),niter,sigma,omega,method);
Ik      = Vk*xk;

% plot

figure;imagesc(reshape(I,n,n),[0 1]);colormap(gray)
axis equal tight; set(gca,'xticklabel',[],'yticklabel',[])

figure;qtplot(Sk);
axis equal tight ij; set(gca,'xticklabel',[],'yticklabel',[]);
set(get(gca,'children'),'linewidth',2)


figure;imagesc(reshape(If,n,n),[0 1]);colormap(gray)
axis equal tight; set(gca,'xticklabel',[],'yticklabel',[])

figure; 
imagesc(reshape(Ik,n,n),[0 1]);colormap(gray);hold on;
axis equal tight;%xlim([20 100]);ylim([30 70]);
set(gca,'xticklabel',[],'yticklabel',[])




%
%savefig(1:4,['../DGCIdoc/figs/' mfilename])