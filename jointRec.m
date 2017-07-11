function [xs,Op] = jointRec(W,p,x0,A,L,lambda,kappa,MaxIter,innerMaxIter,ipStr)
%jointRec Joint reconstruction framework for linear inverse problem
%   The problem Wx = p, is split into two problems by defining
%       x = [1 - h(\phi)]u0 + u1*h(\phi)
% u0 is the background (variational) parameter and u1 is the anomaly 
% parameter which is assumed to be constant. The optimization is done wrt
% to u0 and \phi (a level set function), assuming u1 to be known. The outer
% optimization is wrt to \phi, which is non-linear due to heaviside and
% inner optimization is wrt u0. Hence,
%   min_{u0,\phi} ||Wx - p||^2 = min_{u0} ( min_{\phi} ||Wx - p||^2
%
% Input
%   W       - Linear operator, could be a matrix of size m x n
%   p       - Data, a vector of size m x 1
%   lambda  - regularization parameter
%   kappa   - heaviside width
%   maxit   - maximum number of iterations
%   iPstr   - input structure containing proj_geom and vol_geom
%
% Output
%   xs      - final solution


%% define u0 and u1

global u0 
u1      = 1;                    % anomaly parameter, set to be 1
u0      = zeros(size(W,2),1);   % initialize u0
ng      = sqrt(size(W,2));


%% get RBF Kernel

% Koptions.tau    = 5;        % how coarse the RBF grid should be wrt computational grid
% Koptions.eta    = 4;        % parameter to control the spread of RBF
% Koptions.nouter = 2;        % RBF layers outside compuational domain 
% Koptions.rtype  = 'compact';% RBF type
% Koptions.ltype  = 'L2';     % distance norms for RBF
% 
% x               = 0:(1/(ng-1)):1;    % x-dn vector 
% z               = 0:(1/(ng-1)):1;    % z-dn vector
% 
% [A,nr]          = generateKernel(x,z,Koptions);   
% A               = opMatrix(A);

%% get Regularizer

% n       = [ng ng];  % number of gridpoints
regMode = 1;        % type of regularizer/penalty, 1: Tikhonov for smoothness, 2: penalty for x - x0


%% initialize alpha, level-set parameter

% x0      = -1*ones(nr);
% 
% % make sure to put some positive RBF to create level-set boundary
% x0(floor(nr(1)/2)-3:floor(nr(1)/2),floor(nr(2)/2):floor(nr(2)/2)+3) = 2;    % initial level-set
% 
% x0      = x0(:);

%% minimize over alpha

fh                      = @(x) funcPhi(x,A,u1,W,p,L,regMode,lambda,kappa,innerMaxIter,ipStr);
options                 = optimoptions('fminunc','Algorithm','trust-region','GradObj','on',...
                            'Hessian','user-supplied','HessMult',@(Hinfo,Y)hmfun(Hinfo,Y),'MaxPCGIter',10);
options.Display         = 'iter';
options.MaxIter         = MaxIter;
options.TolFun          = 1e-6;
options.TolX            = 1e-6;
% options.TolPCG          = 1e-2;
% options.PrecondBandWidth= 3;

tic;
[xf,funval,exitflag,Op] = fminunc_new(fh,x0,options);
toc

Op.xf       = xf;
Op.funval   = funval;
Op.exitflag = exitflag;

%% final value
Axf         = A*xf;
hopt.epsi   = 0;
hopt.thr    = kappa*(max(Axf)-min(Axf));
hxf         = heavi(Axf,hopt);
xs          = (1 - hxf).*u0 + hxf*u1;

end


function [f,g,Hinfo] = funcPhi(x,A,u1,W,p,L,regMode,lambda,kappa,innerMaxIter,ipStr)

global u0 fIter

fig         = ipStr.fig;
% [n,~]       = size(A);

%% get h(\phi) and its sensitivity d(\phi)

Ax          = A*x;
hopt.epsi   = 0;    % kappa*(max(Ax) - min(Ax));
[h,~]       = heavi(Ax,hopt);

%% minimize over u_0

% lambda = ((0.9)^fIter) * lambda;
Wu          = [W*opDiag(1-h);sqrt(lambda)*L];

% get data vector 'pu', according to a regularization/penalty
% if regMode == 1
pu      = [p-u1*(W*h);zeros(size(L,1),1)];
% elseif regMode == 2
%     u0est   = 0.25;
%     pu      = [p-u1*(W*h);u0est*ones(n,1)];
% end

% optimize over u0, use lsqlin with bounds constraints
[u0,~]      = lsqr(Wu,pu,1e-3,innerMaxIter,[],[]);



%% define F and data

hopt.epsi   = kappa*(max(Ax) - min(Ax));
[h,d]       = heavi(Ax,hopt);

Um          = opDiag(u1-u0);
D           = opDiag(d);
F           = W*(Um*h);
data        = p - W*u0;


%% function, gradient, Hessian info

f           = 0.5*norm(F - data,2).^2;
% g0          = D*(Um*(W'*(F - data)));
g           = A'*(D*(Um*(W'*(F - data))));
Hinfo       = W*Um*D*A;
% xs          = (1 - h).*u0 + h*u1;

%% plots

% if fig.show
%     fig99 = figure(99);
%     ng = sqrt(size(Ax,1));
%     [zz,xx] = ndgrid(linspace(0,1,ng)); 
%     AxM = reshape(Ax,[ng ng]);
%     % AxD = reshape(d,[ng ng]);
%     % Ztr = (max(Ax) + min(Ax))/2;
%     Ztrm= min(Ax) - (max(Ax)-min(Ax));
%     surfl(xx,zz,AxM); shading interp; hold on;
%     lighting gouraud
%     hfig1 = surf(xx,zz,0*ones(size(AxM)),AxM,'EdgeColor','none');
%     hfig1.FaceAlpha = 0.25;
%     hfig2 = surf(xx,zz,Ztrm*ones(size(AxM)),AxM,'EdgeColor','none');
%     hfig2.FaceAlpha = 0.5;
%     % hfig3 = surf(xx,zz,2*Ztrm*ones(size(AxD)),AxD,'EdgeColor','none');
%     % hfig3.FaceAlpha = 0.5;
%     view(-20,25); set(gca,'GridAlpha',0);set(gca,'color','none');
%     xlabel('x','FontSize',12); ylabel('y','FontSize',12); axis([0 1 0 1 Ztrm max(Ax)]);
%     hold off;
%     pause(0.001);
% %     fig99 = figure(99);
% %     ng = sqrt(size(Ax,1));
% %     subplot(2,2,1); imagesc(reshape(Ax,ng,ng));axis equal tight; axis off;% colorbar; 
% %     title(sprintf('level-set function %.2g',lambda),'Fontsize',10);
% %     subplot(2,2,2); imagesc(reshape(abs(g0),ng,ng));axis equal tight;axis off;% colorbar; 
% %     title(sprintf('level-set sensitivity %.2g',lambda),'Fontsize',10);
% %     subplot(2,2,3); imagesc(reshape(xs,ng,ng));axis equal tight;axis off;% colorbar; 
% %     title(sprintf('reconstructed model %.2g',lambda),'Fontsize',10);
% %     subplot(2,2,4); imagesc(reshape(u0,ng,ng),[0 0.5]);axis equal tight;axis off;% colorbar;
% %     title(sprintf('reconstructed background %.2g',lambda),'Fontsize',10);
% %     pause(0.01);
%     
%     if fig.save
%         fig99.PaperPositionMode = 'auto';
%         saveas(fig99,strcat(fig.path,'evol',num2str(fIter)),'png');
%         fIter = fIter + 1;
%     end
% end
% fIter = fIter + 1;
end


function [Fv] = hmfun(Hinfo,Y)
%hmfun function for Hessian Vector product
% Input
%   Hinfo   - Hessian information from objective function
%   Y       - vector, this comes from optimization procedure
%
% Output
%   Fv      - Matrix-vector product of Hessian and vector

% convert Hinfo matrix to Operator of kind OpSPOT

Fv = (Hinfo)'*(Hinfo*Y);

end