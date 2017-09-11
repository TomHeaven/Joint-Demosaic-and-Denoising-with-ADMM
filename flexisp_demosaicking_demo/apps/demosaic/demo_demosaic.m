% Copyright 2014 NVIDIA Corporation.  All Rights Reserved
% This software is for provided for non-commercial and research purposes only
% Paper: FlexISP: A Flexible Camera Image Processing Framework, SIGGRAPH
% Asia 2014
% WebSite: https://research.nvidia.com/publication/flexisp-flexible-camera-image-processing-framework
% Application: Demosaicking

close all
clear all

addpath('../../3rdparty/BM3D') %%%% download from: http://www.cs.tut.fi/~foi/GCF-BM3D/, if you run into check_order error, you could download the check_order.m from here: https://searchcode.com/codesearch/view/13666439/
addpath('../../3rdparty/LASIP_Image_Restoration_DemoBox_v113'); %%%% download from: http://www.cs.tut.fi/~lasip/2D/
addpath('../../3rdparty/Gabriel_Peyre'); %%%% copy grad.m, getoption.m, div.m to this folder. download from: https://www.ceremade.dauphine.fr/~peyre/codes/
addpath('../../3rdparty/Matlab'); %%%%% if you have Partial Differential Equation/Signal Processing toolbox, you are good. Otherwise put dst.m, dct.m and kaiser.m to this folder. 
addpath('../../utils')
addpath('../../core')
addpath('../../core/priors')
addpath('../../data/demosaic/')




%% Configuration
% cfa pattern
conf_cfa = 'rggb';

% configuration for primal dual optimization framework
conf_pdopt_sigma_scale = 12.5;
conf_pdopt_lambda_reg = 1.0;                  
conf_pdopt_max_iters = 30;

conf_lambda_channels = [[2000,   [0.0,  0.0, 0.0]]; ... %R
                        [2000,   [0.0,  0.0, 0.0]]; ... %G
                        [2000,   [0.0,  0.0, 0.0]]];    %B

%% Input images
% read in image
I = im2double(imread('input.tif'));

% generate raw format
I_raw_input = generate_bayer(I, conf_cfa);

%% Prepare data for optimizer format

% observed signal
observed_input(:,:,1) = I_raw_input;
% mask for the observed signal.
observed_input(:,:,2) = ones(size(I_raw_input));
% initial guess
initial_guess = double(demosaic( uint8(observed_input(:,:,1) * 255.0), conf_cfa)) / 255.0;

% generate mask for the input data (for backprojection purpose)
[red_mask, green_mask, blue_mask] = generate_bayer_mask(size(I(:, :, 1)), conf_cfa);
channel_mask = cat(3, red_mask, green_mask, blue_mask);

% blur kernel
blur_size = 5; %has to be an odd number !
kernel_blur = cell(0); % blur kernel for each channel
channel_patch = cell(0);
for ch = 1:size(I,3)
    % no blur
    kernel_blur{ch} = zeros([blur_size, blur_size]);
    kernel_blur{ch}( floor(blur_size/2) + 1, floor(blur_size/2) + 1 ) = 1; 
    
    % channel patch
    channel_patch(ch).Image = observed_input(:, :, 1) .* channel_mask(:, :, ch);
    channel_patch(ch).K = kernel_blur{ch};
    channel_patch(ch).sat_mask = ones(size(observed_input));
end

% initial guess
initial_channel_patch = channel_patch;
for ch = 1:size(I,3)
    initial_channel_patch(ch).Image = initial_guess(:, :, ch);
end

%% Run optimization
fprintf('data ready! run optimization\n');
result = pd(channel_patch, initial_channel_patch, ...
            conf_lambda_channels, conf_pdopt_lambda_reg, conf_pdopt_sigma_scale, conf_pdopt_max_iters , ...
            1e-6, '.', 'brief', I);

result_demosaic = cat(3, result(1).Image, result(2).Image, result(3).Image);

%% compute psnr
% remove the artifact border
border = 5;
sr = border; er = size(I, 1) - border;
sc = border; ec = size(I, 2) - border;

% red channel
psnr_r = comppsnr(result_demosaic(sr:er, sc:ec, 1), I(sr:er, sc:ec, 1));
% green channel
psnr_g = comppsnr(result_demosaic(sr:er, sc:ec, 2), I(sr:er, sc:ec, 2));
% blue channel
psnr_b = comppsnr(result_demosaic(sr:er, sc:ec, 3), I(sr:er, sc:ec, 3));

fprintf('PSNR = [%.2f %.2f %.2f]\n', psnr_r, psnr_g, psnr_b);

figure(1);
imshow(result_demosaic);
pause(0.1)

return