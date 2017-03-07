function f_seg = segment(f_org,grey_values)
%--------------------------------------------------------------------------
% f_thres = segment(f_org,grey_values)
%               
% Computes the segmentation of an image by rounding grey values to the
% nearest segmentation value.
%
% Input:
% f_org         The non-segmented image;
% grey_values   The segmentation values.
%
% Output:
% f_seg     	The segmented image.
%--------------------------------------------------------------------------

thresholds = grey_values(1:end-1) + diff(grey_values)/2;

f_seg = ones(size(f_org))*grey_values(1);

for n = 2:length(grey_values)

    f_seg(thresholds(n-1) < f_org) = grey_values(n);

end
    
end