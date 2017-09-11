% Copyright 2014 NVIDIA Corporation.  All Rights Reserved
% This software is for provided for non-commercial and research purposes only
% Paper: FlexISP: A Flexible Camera Image Processing Framework, SIGGRAPH
% Asia 2014
% WebSite: https://research.nvidia.com/publication/flexisp-flexible-camera-image-processing-framework

function [ psnr, mse ] = comppsnr( x, y )
%PSNR Summary of this function goes here
%   Detailed explanation goes here

xl = x(~isnan(x) & ~isnan(y));
yl = y(~isnan(x) & ~isnan(y));

maxval = 1.0;

if (max(xl(:)) > 10.0) && (max(yl(:)) > 10.0)
    maxval = 255;
elseif max(xl(:)) > 10.0 
    xl = xl / 255;
elseif max(yl(:)) > 10.0
    yl = yl / 255;
end

diff2 = (xl-yl).^2;
mse = mean(diff2(:));
psnr = 10.0 * log10(maxval*maxval / mse);
end
