function edges = boundary(im_org)
%--------------------------------------------------------------------------
% edges = boundary(im_org)
%
% Determine the edges of a segmented image according to a 8-conntected
% neigbhbourhood.
%
% Input:
% im_org        The segmented image.
%
% Output:
% edges         An image which is zero if the pixel is not an edge and
%               greater than zero if it is, the value indicates how many
%               neighbouring pixels have a different grey value.
%--------------------------------------------------------------------------

if ndims(im_org) == 2
    edges = boundary2D(im_org);
elseif ndims(im_org) == 3
    edges = boundary3D(im_org);
else
    error('The object has invalid dimensions!');
end

end

function edges = boundary2D(im_org)
kernel = [1 1 1; 1 1 1; 1 1 1];


im_large = padarray(im_org, [1,1]);
% im_large(2:end-1, 2:end-1) = im_org;

edges = zeros(size(im_org));

for s = -1:1
    
    for t = -1:1
        
        if kernel(s+2, t+2) ~= 0
           
            Temp = abs(im_large(2:end-1, 2:end-1) - im_large(2+s:end-1+s, 2+t:end-1+t));
            edges(Temp > eps) = edges(Temp > eps) + 1;
        
        end
        
    end
    
end

end


function edges = boundary3D(S)
			vol_geom = astra_create_vol_geom(size(S,2),size(S,1),size(S,3));
			data_id = astra_mex_data3d('create', '-vol', vol_geom, S);

			cfg = astra_struct('DARTMASK3D_CUDA');
			cfg.SegmentationDataId = data_id;
			cfg.MaskDataId = data_id;
			cfg.option.GPUindex = 0;
% 			cfg.option.Connectivity = 6;
            cfg.option.Connectivity = 26;
			cfg.option.Threshold = 1;
			cfg.option.Radius = 1;
			
			alg_id = astra_mex_algorithm('create',cfg);	
			astra_mex_algorithm('iterate',alg_id,1);
			edges = astra_mex_data3d('get', data_id);
		
			astra_mex_algorithm('delete', alg_id);
			astra_mex_data3d('delete', data_id);
			
			% random
% 			RandomField = double(rand(size(S)) < this.random);

			% combine
% 			Mask = or(Mask, RandomField);
% objectPadded = padarray(object, [1,1,1]);
% 
% edges = zeros(size(object));
% 
% for x = -1:1
%     for y = -1:1
%         for z = -1:1
%             idx = abs(objectPadded((2:end-1)+x, (2:end-1)+y, ...
%                 (2:end-1)+z) - object) > eps;
%             edges(idx) = edges(idx) + 1;
%         end
%     end
% end

end
