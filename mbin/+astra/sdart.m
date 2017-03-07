function [recon_seg,recon_grey,times,pixelError] = sdart(W,projections,volSize,greyValues,method,initial_arm_it,arm_it,dart_it,mode, lambda, x0, x_true, nullspace)
%SDART Soft Discrete Algebraic Reconstruction Technique
%   RECON_SEG = SDART(W,PROJECTIONS,VOLSIZE,GREY_VALUES) computs the
%   SDART reconstruction from the PROJECTIONS at an image size of VOLSIZE.
%   This method tries to solve the system:
%       [ W ]      [ p  ]
%       [ D ]  x = [ Dv ]
%   where D is a regularization matrix and v is the segmentation of the
%   current reconstrucion. W is the projection operator and should be a
%   matrix or Spot. The PROJECTIONS can be a vector (consistent with W) or
%   have the dimensions of a sinogram or a stack of image. The 
%   reconstructed image is constrained to the provided GREYVALUES.
%   The output RECON_SEG reconstruction has been segmented.
%
%   RECON_SEG = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,METHOD) also
%   pass the algebraic reconstruction method (ARM) as a string, if METHOD
%   is empty the default is used, 'cgls'. Other possibilities for METHOD 
%   are:
%     'cgls', 'cgne', 'lsqr', 'lsmr', 'scaled_cgls', 'scaled_cgne',
%     'scaled_lsqr', 'scaled_lsmr'.
%
%   RECON_SEG = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,METHOD, 
%   INITIAL_ARM_IT) also defines the number of iterations used for the 
%   initial ARM reconstruction, if INITIAL_ARM_IT is empty the default is 
%   used, 50.
%
%   RECON_SEG = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,METHOD,
%   INITIAL_ARM_IT,ARM_IT) also pass the number of intermediate ARM
%   iterations. If ARM_IT is empty, the default will be used, 20.
%
%   RECON_SEG = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,METHOD,
%   INITIAL_ARM_IT,ARM_IT,DART_IT) also pass the number of SDART
%   iterations. If DART_IT is empty, the default will be used, 20.
%
%   RECON_SEG = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,METHOD,
%   INITIAL_ARM_IT,ARM_IT,DART_IT,MODE) als pass the regularization type
%   MODE for matrix D. If MODE is empty the default will be used, 3.
%   Possible values for MODE are:
%       1 - Mimics classic DART
%       2 - Same as 1, but "free" pixels are slightly restricted
%       3 - Neighbor criterion, pixels with many similar neighbors are more
%           restricted than those having many different neighbors
%       4 - D_ii depends on the total variation of the neighbors of pixel i
%       5 - D_ii is based on the distance from pixel i to the nearest pixel
%           with a different grey value
%       6 - D_ii is based on the difference between the pixel value and its
%           segmented value.
%       7 - D_ii is based on the neighbor similarity of pixels within a
%           disk centered around pixel i. The radius of this disk shrinks
%           with the iteration number

%
%   RECON_SEG = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,METHOD,
%   INITIAL_ARM_IT,ARM_IT,DART_IT,MODE,LAMBDA) also pass the regularization
%   parameter LAMBDA that influences the regularization matrix D. If LAMBDA
%   is not chosen right, the method might not converge!
%
%   RECON_SEG = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,METHOD,
%   INITIAL_ARM_IT,ARM_IT,DART_IT,MODE,LAMBDA,X0) also pass the initial
%   solution X0. The default is a blank image (all zeros).
%
%   RECON_SEG = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,METHOD,
%   INITIAL_ARM_IT,ARM_IT,DART_IT,MODE,LAMBDA,X0,X_TRUE) also pass the 
%   ground truth image X_TRUE. This is used to compute the pixel error 
%   which can be useful in a phantom study.
%
%   RECON_SEG = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,METHOD,
%   INITIAL_ARM_IT,ARM_IT,DART_IT,MODE,X0,X_TRUE,NULLSPACE) determine if 
%   the intermediate ARM update starts with a blank image (all zeros, 
%   NULLSPACE==TRUE) or if the previous grey value reconstruction is used
%   to initate the iterations (NULLSPACE==FALSE, default).
%
%   [RECON_SEG,RECON_GREY] = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,...)
%   also returns the reconstruction RECON_GREY that is not segmented.
%
%   [RECON_SEG,RECON_GREY,TIMES] = SDART(W,PROJECTIONS,VOLSIZE,GREYVALUES,
%   ...) also returns the wallclock times TIMES that have passed in each
%   iteration.
%
%   [RECON_SEG,RECON_GREY,TIMES,PIXELERROR] = SDART(W,PROJECTIONS,VOLSIZE,
%   GREYVALUES,...) also returns the PIXELERRORs of the segmented
%   reconstruction at each iteration. Note that the ground truth should be
%   provided.
%
%   See also DART, ADART.

%   Copyright 2012, 2013 Frank Tabak, and Folkert Bleichrodt

% input parsing
if nargin < 4  ,  error('Not enough input arguments!');         end
if nargin < 5  || isempty(method)        , method = 'cgls';     end
if nargin < 6  || isempty(initial_arm_it), initial_arm_it = 50; end
if nargin < 7  || isempty(arm_it)        , arm_it  = 20;        end
if nargin < 8  || isempty(dart_it)       , dart_it = 20;        end
if nargin < 9  || isempty(mode)          , mode = 3;            end
if nargin < 10 || isempty(lambda)        , lambda = 1;          end
if nargin < 11 || isempty(x0)            , x0 = [];             end
if nargin < 12 || isempty(x_true)        , x_true = [];         end
if nargin < 13 || isempty(nullspace)     , nullspace = false;   end

groundtruthAvailable = ~isempty(x_true);

N = size(W,2);

if ~iscolumn(projections)
    warning('Assuming column-wise storage for sinogram!');
end

p = projections(:);

times = zeros(dart_it,1);

% Initial ARM Reconstruction
display('Initial Reconstruction');

if isempty(x0)
    tic;
    recon_grey = astra.applyArm(method, W, p, initial_arm_it);
    toc;
else
    recon_grey = x0(:);
end

fprintf('iter\t\trel. residual');

if groundtruthAvailable
    fprintf('\tpixel error');
    pixelError = zeros(dart_it,1);
end

fprintf('\n--------------------------------------------\n');

t_start = tic;

% Start DART iterations
for i = 1:dart_it-1

    % Segment continuous image
    recon_seg = astra.utils.segment(recon_grey,greyValues);

    fprintf('%d of %d,\t %f', i, dart_it, norm(W*recon_seg - p)/norm(p));
    if groundtruthAvailable
        pixelError(i) = sum(recon_seg(:) ~= x_true(:));
        fprintf('\t%d', pixelError(i));
    end
    fprintf('\n');
    
    imshow(reshape(recon_seg,volSize), []);
    drawnow;
    
    % Identify boundary and free pixels
    free  = astra.utils.boundary(reshape(recon_seg,volSize));
    fixed = ~free;
    free_columns  = find(free > 0);
    
    % Create new (regularized) system
    switch (mode)
        case 1
            % Fixed variables get very high d_j and v is segmented solution as
            % to steer the system to a solution with fixed variable equal to the
            % segmented value. Free variables have d_j = 0
            D     = lambda * 1e6 * diag(sparse(reshape(fixed,N,1)));
            W_reg = [W; D];
            v     = recon_seg;
            v(free_columns) = 0;
            p_reg = [p; D*v];
            
        case 2
            % Fixed variables get very high d_j and v is segmented solution as
            % to steer the system to asolution with fixed variable equal to the
            % segmented value. Free variables have d_j = 1. Thus there is a
            % Tikhonov regularization on the free pixels.
            d_i   = 10^6*ones(N,1);
            d_i(free_columns) = 1;
            D     = diag(sparse(d_i));
            W_reg = [W; D];
            v     = recon_seg;
            v(free_columns) = 0;
            p_reg = [p; D*v];
            
        case 3
            % The diagonal elements depend on the number of neighours a 
            % pixel has.
            D      = diag(sparse(100*(1./3.^(free(:)))));
            D      = lambda*D;

            W_reg = [W; D];
            v     = recon_seg;
            p_reg = [p; D*v];
            
        case 4
            % The value of d_j is now depent upon the total variation of the
            % grey values (not segmented!) of the neighbours of pixel i.
            im_grey     = reshape(recon_grey,volSize);
            im_grey_pad = padarray(im_grey,[1 1],'replicate');
            tot_var_big = abs(conv2(im_grey_pad,[1 1 1; 1 -8 1; 1 1 1],'same'));
            tot_var     = tot_var_big(2:end-1,2:end-1);

            D     = lambda * diag(reshape(sparse(1e2./(10.^(tot_var))),N,1));
            W_reg = [W; D];
            v     = recon_seg;
            p_reg = [p; D*v];
            
        case 5
            % Pixels are weighted based on their distance to the nearest
            % boundary of two homogeneous regions.
            distToBoundary = double(bwdistMulti(reshape(recon_seg, volSize))).^2;
            D     = lambda * diag(sparse(distToBoundary(:)));
            W_reg = [W; D];
            v     = recon_seg;
            p_reg = [p; D*v];
            
        case 6
            W_reg  = [W; lambda*speye(N)];
            v      = recon_seg;
            p_reg  = [p; lambda*speye(N)*v];
            
        case 7
            nnb   = computeNeighborVariety(reshape(recon_seg,volSize), max(10-i,1));
            nnb   = nnb/max(nnb(:));
            % rescale
            nnb = 100*3.^(-8*(nnb));

            D     = diag(sparse(nnb(:)));
            W_reg = [W; D];
            v     = recon_seg;
            p_reg = [p; D*v];
            
        otherwise
            error('I did not recognize that mode for soft DART?');          
    end

    % Apply ARM to the new system
    if nullspace
        recon_grey = astra.applyArm(method, W_reg, p_reg, arm_it);
    else
        recon_grey = astra.applyArm(method, W_reg, p_reg, arm_it, recon_grey);
    end
    
    times(i) = toc(t_start);
end

% Segement final reconstruction
recon_seg = astra.utils.segment(recon_grey, greyValues);

fprintf('%d of %d,\t %f', dart_it, dart_it, norm(W*recon_seg - p)/norm(p));
if groundtruthAvailable
    pixelError(dart_it) = sum(recon_seg(:) ~= x_true(:));
    fprintf('\t%.2e', pixelError(dart_it));
end
fprintf('\n');

times(dart_it) = toc(t_start);

end % sdart

%% Extra utility functions
function y = bwdistMulti(x)
% BWDISTMULTI
%   Y = BWDISTMULTI(X) for the image X, for each pixel, the Euclidean
%   distance to the nearest nonzero pixel is computed.

% assume the image has been segmented
values = unique(x);
nVal   = numel(values);

y = zeros(size(x));

for iVal = 1:nVal
    y = y + bwdist(x ~= values(iVal));
end

end % bwdistMulti


function y = computeNeighborVariety(img, radius)
%COMPUTENEIGHBORVARIETY
%   Y = COMPUTENEIGHBORVARIETY(IMG,RADIUS) counts for each pixel in
%   image IMG the number of neighbors inside a disk with radius RADIUS,
%   that have a grey value that differs from the center pixel.

if ismatrix(img)
    y = computeNeighborVariety2D(img, radius);
elseif ndims(img == 3)
    y = computeNeighborVariety3D(img, radius);
else
    error('Dimensions of the object are invalid!');
end

end % computeNeigborVariety


function y = computeNeighborVariety2D(img, radius)
[mRows, nCols] = size(img);

% intialize y
y = zeros(size(img));

% pad image
img = padarray(img, [radius,radius], 0);

% construct "convolution" kernel
kernel = circConv(radius);
[mKernel, nKernel] = size(kernel);

for iRow = 1:mRows
    for jCol = 1:nCols
        greyVal       = img(iRow+radius, jCol+radius);
        subImage      = img(iRow:mKernel+iRow-1, jCol:nKernel+jCol-1);
        % Count neighbors that have different grey value
        y(iRow, jCol) = sum(subImage(kernel) ~= greyVal);
    end
end

end % computeNeigborVariety2D


function y = computeNeighborVariety3D(img, radius)

% initialize y
y = zeros(size(img));

[mRow, nCol, pSlice] = size(img);

% pad image
img = padarray(img, [radius,radius,radius]);

% construct "convolution" kernel
kernel = sphereConv(radius);
[mKernel, nKernel, pKernel] = size(kernel);

for iRow = 1:mRow
    for jCol = 1:nCol
        for kSlice = 1:pSlice
            greyVal  = img(iRow+radius, jCol+radius, kSlice+radius);
            subImage = img(iRow:mKernel+iRow-1, jCol:nKernel+jCol-1, ...
                kSlice:pKernel+kSlice-1);
            % count neighbors that have different grey value
            y(iRow, jCol, kSlice) = sum(subImage(kernel) ~= greyVal);
        end
    end
end

end % computeNeighborVariety3D


function kernel = circConv(radius)
%CIRCCONV Creates a circular convolution kernel
%   KERNEL = CIRCCONV(RADIUS) generates a (RADIUS+1) x (RADIUS+1) matrix
%   with elements 0 and 1. The nonzeros form a circle with radius RADIUS

[x,y] = meshgrid(-radius:radius, -radius:radius);

kernel = (x.^2+y.^2) <= radius*radius;

end % circConv


function kernel = sphereConv(radius)

[x,y,z] = meshgrid(-radius:radius, -radius:radius, -radius:radius);

kernel = (x.^2+y.^2+z.^2) <= radius*radius;

end % sphereConv
