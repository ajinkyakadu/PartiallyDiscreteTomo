function [reconstruction, segmentation, pixelError, times, allrecons] = dart(projections, proj_geom, vol_geom, greyValues, init_arm_it, arm_it, dart_it, method, fix_prob, x0, im_true)
%function [reconstruction, segmentation, pixelError, times, allrecons] = dart(projections, proj_geom, vol_geom, greyValues, init_arm_it, arm_it, dart_it, method, fix_prob, x0, im_true)

if (nargin < 4) , error('Not enough input arguments!');           end
if (nargin < 5  || isempty(init_arm_it)), init_arm_it = 50;       end
if (nargin < 6  || isempty(arm_it))     , arm_it   = 20;          end
if (nargin < 7  || isempty(dart_it))    , dart_it  = 20;          end
if (nargin < 8  || isempty(method))     , method   = 'SIRT_CUDA'; end
if (nargin < 9  || isempty(fix_prob))   , fix_prob = 1;           end
if (nargin < 10 || isempty(x0))         , x0 = [];                end
if (nargin < 11 || isempty(im_true))    , im_true = [];           end

if nargout >= 5
    allrecons = zeros(vol_geom.GridRowCount * vol_geom.GridColCount,dart_it);
end

groundtruthAvailable = ~isempty(im_true);

if groundtruthAvailable
    pixelError = zeros(dart_it,1);
end

times = zeros(dart_it,1);

if iscolumn(projections)
    projections = reshape(projections,[numel(proj_geom.ProjectionAngles),proj_geom.DetectorCount]);
end

is2D = ~isfield(vol_geom, 'GridSliceCount');

if ~isempty(x0)
    error('Cannot use custom initial solution.');
end

D = DARTalgorithm(projections, proj_geom);
D.t0 = init_arm_it;
D.t = arm_it;

if is2D
    D.tomography = IterativeTomography();
else
    D.tomography = IterativeTomography3D();
end

D.tomography.vol_geom     = vol_geom;
D.tomography.method       = method;
D.tomography.inner_circle = 'no';
D.tomography.proj_type    = 'strip';


% compute tresholds
thresHolds = greyValues(1:end-1) + diff(greyValues)/2;

D.segmentation.rho     = greyValues;
D.segmentation.tau     = thresHolds;

D.smoothing.b          = 0.50;
D.smoothing.gpu        = 'yes';
% D.smoothing.full3d     = 'yes';
% D.smoothing.gpu_core   = gpu_core;
 
D.masking.random       = 1 - fix_prob;

if is2D
    D.masking.conn = 8;
else
    % For 3D 26 is also an option
%     D.masking.conn = 26;   
   D.masking.conn = 6;
end

% D.output.directory     = 'outDart/';
% D.output.pre           = 'data';
% D.output.save_images   = 'no';
% D.output.save_results  = {'stats', 'settings', 'S', 'V'};
% D.output.save_interval = 1;
% D.output.verbose       = 'no';

D.statistics.proj_diff = 'no';

D.initialize();

disp([D.output.directory D.output.pre]);

% global im

t_start = tic;
% for i = 1:dart_it
%     D.iterate(1);
%     if nargin == 11
%         pixelError(i) = sum(sum(D.S ~= im_true));
%         fprintf('%d\n', pixelError(i));
%     end
%     if is2D
%         imshow(D.S,[]);
% %         imwrite(D.S(257:end-256,257:end-256), sprintf('../img/anim_bone/segmentation/frame%03d.png', i));
% %         err = im - D.S;
%         % scaling
% %         err_min = min(greyValues) - max(greyValues);
% %         err_max = max(greyValues) - min(greyValues) - err_min;
% %         err = (err - err_min)/err_max;
% %         min(err(:))
% %         max(err(:))
% %         imshow(err(257:end-256,257:end-256), [min(greyValues) - max(greyValues), max(greyValues) - min(greyValues)]);
% %         
% %         err_min = min(greyValues) - max(greyValues);
% %         err_max = max(greyValues) - min(greyValues);
% %         
% %         err(err <= err_min) = err_min;
% %         err(err >= err_max) = err_max;
% %         
% %         show(err(257:end-256,257:end-256));
% %         
% %         imwritesc(err(257:end-256,257:end-256), sprintf('../img/anim_bone/error/frame%03d.png', i));
%         drawnow;
%         if nargout >= 5
%             allrecons(:,i) = column(D.S);
%         end
%     else
%         % show middle slice
%         imshow(D.S(:,:,round(size(D.S,3)/2))',[]);
%         drawnow
%     end
%     times(i) = toc(t_start);   
% end
    
D.iterate(dart_it);

segmentation = D.S;

reconstruction = D.V;
