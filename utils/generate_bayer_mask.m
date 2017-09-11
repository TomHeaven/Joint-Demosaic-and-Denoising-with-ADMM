% Copyright 2014 NVIDIA Corporation.  All Rights Reserved
% This software is for provided for non-commercial and research purposes only
% Paper: FlexISP: A Flexible Camera Image Processing Framework, SIGGRAPH
% Asia 2014
% WebSite: https://research.nvidia.com/publication/flexisp-flexible-camera-image-processing-framework

function [ r_mask, g_mask, b_mask ] = generate_bayer_mask( sz, pattern )
% generate color mask for bayer

r_mask = zeros(sz);
g_mask = zeros(sz);
b_mask = zeros(sz);

if strcmp(pattern, 'rggb')
    % green
    g_mask(1:2:end, 2:2:end) = 1;
    g_mask(2:2:end, 1:2:end) = 1;

    % red
    r_mask(1:2:end, 1:2:end) = 1;

    % blue
    b_mask(2:2:end, 2:2:end) = 1;
elseif strcmp(pattern, 'grbg')
    % green
    g_mask(1:2:end, 1:2:end) = 1;
    g_mask(2:2:end, 2:2:end) = 1;

    % red
    r_mask(1:2:end, 2:2:end) = 1;

    % blue
    b_mask(2:2:end, 1:2:end) = 1;
end

end

