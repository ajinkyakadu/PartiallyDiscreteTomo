function [fig] = plotRec(x,W,p,im,diff)

reconstruction = reshape(x, W.vol_size);

fig = figure();
subplot(2,2,1);
imagesc(reconstruction,[0 1]); axis equal tight; axis off;
title('Reconstruction');

subplot(2,2,2);
imagesc(im,[0 1]); axis equal tight; axis off;
title('Ground truth');

% The transpose of the operator corresponds to the backprojection.
backProjection = W'*p;
subplot(2,2,3);
imagesc(reshape(backProjection, W.vol_size)); axis equal tight; axis off;
title('Backprojection');

subplot(2,2,4);
imagesc(reshape(diff,W.vol_size),[-1 1]); axis equal tight; axis off;
title('final difference');

end
