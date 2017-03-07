function x = chambollePock(A, TV, b, maxit, lambda, nonnegative, x0, visualize)
%CHAMBOLLEPOCK TV-minimization with Chambolle-Pock
%   X = CHAMBOLLEPOCK(A,TV,B) aims at solving:
%      minimize_x ||A*x-B|| + \lambda* ||TV(x)||_1
%   Either A and/or TV can be function handles that operate as:
%      A(x,1)  -  returns A*x
%      A(x,2)  -  returns A'*x
%   As an example one can use the following matrix for the TV operator:
%      n = size(A,2);
%      D = spdiags([[-ones(n-1,1);0],ones(n,1)], [0,1], n, n);
%      TV = [kron(speye(n), D); kron(D, speye(n))];
%
%   X = CHAMBOLLEPOCK(A,TV,B,MAXIT) stops the algorithm when MAXIT
%   iterations are reached, default MAXIT = 20.
%
%   X = CHAMBOLLEPOCK(A,TV,B,MAXIT,LAMBDA) pass the TV weight LAMBDA. The
%   default is LAMBDA = 10;
%
%   X = CHAMBOLLEPOCK(A,TV,B,MAXIT,LAMBDA,NONNEGATIVE) specifies if
%   nonnegativity constraints should be applied to the solution.
%
%   X = CHAMBOLLEPOCK(A,TV,B,MAXIT,LAMBDA,NONNEGATIVE,X0) pass initial
%   guess X0.
%
%   X = CHAMBOLLEPOCK(A,TV,B,MAXIT,LAMBDA,NONNEGATIVE,X0,VISUALIZE)
%   determines the level of detail in output. If VISUALIZE is 0, no output
%   is given. If it is 1, error estimates are presented. If it is larger,
%   an image of the solution is shown.

%   Copyright 2013, Folkert Bleichrodt (CWI, Amsterdam)
%
%   This code is derived from the article:
%   [1] Emil Y. Sidky, Jakob H. Jorgensen and Xiaochuan Pan, "Convex
%   optimization problem prototyping for image reconstruction in computed
%   tomography with the Chambolle-Pock algorithm", 2012
%   Download from: \url{http://arxiv.org/abs/1111.5632}

if nargin < 4 || isempty(maxit)      , maxit  = 20;         end
if nargin < 5 || isempty(lambda)     , lambda = 10;         end
if nargin < 6 || isempty(nonnegative), nonnegative = false; end
if nargin < 7 || isempty(x0)         , x0 = 0*A'*b;         end
if nargin < 8 || isempty(visualize)  , visualize   = 1;     end

% Check for explicit matrices
% Note it is slightly dangerours to assume that "is not a function_handle"
% implies "is numeric", however, we do this to support Spot, without doing
% so explicitly ;).
explicitA  = ~isa(A,  'function_handle');
explicitTV = ~isa(TV, 'function_handle');

% TODO: implement function handles!
% if explicitA
%     x0 = 0*A'*b;
% else
%     x0 = 0*A(b,2)';
% end


% inital vectors
u = x0;
p = zeros(size(b,1),1);
q = zeros(size(TV,1),1);

n = numel(x0);
dim = size(TV,1) / size(TV,2);

L = tvmin_powerMethod(A, TV);

% L1 = powerMethod2(A)
% L2 = powerMethod2(TV)

% parameters
tau = 1/L;
sigma = 1/L;
theta = 1;

if visualize > 0
    fprintf('Iter\tresiudal\tprimal-dual gap\n');
    fprintf('---------------------------------------\n');
end

u_bar = u;

for i = 1:maxit
    p = (p + sigma*(A*u_bar - b))/(1+sigma);
    q = q + sigma*TV*u_bar;
  
    q = (lambda*q)./max(lambda,abs(q));
    
    % TODO: make more memory efficient
    % with min constraint
    if nonnegative
        u_new = max(u - tau*A'*p - tau*TV'*q,0);
    else
        u_new = u - tau*A'*p - tau*TV'*q;
    end

    u_bar = u_new + theta*(u_new-u);
    u = u_new;
    
    if visualize > 0
        res = norm(A*u - b) + lambda*norm(TV*u,1);
        gap = norm(A*u - b) + lambda*norm(TV*u,1) + .5*norm(p)^2 + p'*b;
        fprintf('%d\t%e\t%e\n', i, res, gap);
    end    
    if visualize > 1
        if dim == 2
            show(u);
        elseif dim == 3
            nn = round(n^(1/3));
            if nn^3 == n
                U = reshape(u, [nn,nn,nn]);
                show(U(:,:,round(nn/2)));
            end
        end
    end

end

x = u;

end % chambollePock


function L = tvmin_powerMethod(W, D)
%TVMIN_POWERMETHOD Power method for TV-min tomo system matrix
%   L = TVMIN_POWERMETHOD(W, D) computes the largest singular value L of 
%   the matrix [W,D].

niter = 5;

% random initial image
x = rand(size(W,2),1);
y = W*x;
z = D*x;

for i = 1:niter
    % power iteration
    x = W'*y + D'*z;
    % normalize
    x = x./norm(x);
    % One less matvec per iteration,
    % but double the memory usage.
    y = W*x;
    z = D*x;
    L = sqrt(y'*y + z'*z);
end

end

