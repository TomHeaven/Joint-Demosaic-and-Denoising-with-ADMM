% Copyright 2014 NVIDIA Corporation.  All Rights Reserved
% This software is for provided for non-commercial and research purposes only
% Paper: FlexISP: A Flexible Camera Image Processing Framework, SIGGRAPH
% Asia 2014
% WebSite: https://research.nvidia.com/publication/flexisp-flexible-camera-image-processing-framework

function [ raw ] = generate_bayer( I, pattern)
% generate synthetic CFA format with a given image and pattern.

%% check if I is normalized
if((max(I(:)) > 1.0) && (min(I(:)) < 0.0))
    error('the input image is not normalized!');
end

%% generate mask for each channel
[red_mask, green_mask, blue_mask] = generate_bayer_mask(size(I(:, :, 1)), pattern);

%% generate final rggb format
raw = I(:, :, 1) .* red_mask + I(:, :, 2) .* green_mask + I(:, :, 3) .* blue_mask;

end

