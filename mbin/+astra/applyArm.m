function [f_grey] = applyArm(arm, W, p, nIter, x0)
%APPLYARM Applies an Algebraic Reconstruction Method
%   F_GREY = APPLYARM(ARM, W, P) tries to solve the system W*X = P for the
%   unkowns X. The Algebraic Reconstruction Method as passed in the string
%   ARM is used. The resulting reconstrucution is passed as an image
%   F_GREY.
%
%   F_GREY = APPLYARM(ARM, W, P, NITER) also pass the number of iterations
%   NITER after which the algorithm will be stopped. If NITER is [] then
%   APPLYARM uses the default, 100.
%
%   F_GREY = APPLYARM(ARM, W, P, NITER, X0) also passes the initial image
%   X0. If X0 is [] then APPLYARM uses the default, zeros(N,1), with N is
%   size(W,2);

if nargin < 3, error('Not enoug input arguments!'); end

N = size(W,2);

if nargin < 4 || isempty(nIter), nIter = 100;       end
if nargin < 5 || isempty(x0)   , x0 = zeros(N,1);   end


switch arm
    case 'sirt'
        f_grey = astra.arm.sirt(W, x0, p, nIter);
    case 'cgls'
        f_grey = astra.arm.cgls(W,p,[],[],nIter,[],x0);
    case 'cgne'
        f_grey = astra.arm.cgne(W, p, x0, nIter);
    case 'lsqr'
        [f_grey, flag] = lsqr(W, p, [], nIter, [], [], x0);
    case 'lsmr'
        f_grey = astra.arm.lsmr(W, p, [], [], [], [], nIter);
    case 'scaled_cgls'
        f_grey = astra.arm.scaled_cgls(W, p, x0, nIter);
    case 'scaled_cgne'
        C = diag(sparse(1./sum(W)));
        R = diag(sparse(1./sum(W,2)));
        C(isinf(C)) = 0;
        R(isinf(R)) = 0;
        temp   = astra.arm.cgne(sqrt(R)*(W*sqrt(C)), sqrt(R)*p, x0, nIter);
        f_grey = sqrt(C)*temp;
    case 'scaled_lsqr'
        C = diag(sparse(1./sum(W)));
        R = diag(sparse(1./sum(W,2)));
        C(isinf(C)) = 0;
        R(isinf(R)) = 0;
        temp   = astra.arm.lsqrSOL(size(sqrt(R)*(W*sqrt(C)),1), size(sqrt(R)*(W*sqrt(C)),2), sqrt(R)*(W*sqrt(C)), sqrt(R)*p, [], [], [], [], nIter, []);
        f_grey = sqrt(C)*temp;
    case 'scaled_lsmr'
        C = diag(sparse(1./sum(W)));
        R = diag(sparse(1./sum(W,2)));
        C(isinf(C)) = 0;
        R(isinf(R)) = 0;
        temp   = astra.arm.lsmr(sqrt(R)*(W*sqrt(C)), sqrt(R)*p, [], [], [], [], nIter, []);
        f_grey = sqrt(C)*temp;
    case 'cgnr'
        f_grey = astra.arm.cgnr(W, p, nIter, x0);
    otherwise
        error('Please pass a valid arm name!');
end

end

