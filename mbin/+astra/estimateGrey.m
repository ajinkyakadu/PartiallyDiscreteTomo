function greyValue = estimateGrey(proj_geom, vol_geom, projections, USP, range)

projections = projections(:);

W = astra.opTomo('cuda', proj_geom, vol_geom);

nRange = numel(range);

errorVec = zeros(nRange,1);

for iRange = 1:nRange
    % grey value
    g = range(iRange);
    
    projG = projections - W*(g*USP(:));
    
    % reconstruct the complement of USP
    [recon_id, recon] = astra_create_reconstruction_cuda('SIRT_CUDA', proj_geom, vol_geom, reshape(projG,W.proj_size), 100, 'yes', double(~USP), 'no', 0, 'no', 255);
    astra_mex_data2d('delete', recon_id);
    show(recon + g*USP);
    errorVec(iRange) = sqrt(sumsqr(projG - W*recon(:)));
end

[~, i] = min(errorVec);
greyValue = range(i);

if (i == 1)
    warning('Minimum might be smaller than the provided range!');
end

if (i == numel(errorVec))
    warning('Minimum might be larger than the provided range!');
end

end % estimateGrey
