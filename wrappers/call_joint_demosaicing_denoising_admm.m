function RGB = call_joint_demosaicing_denoising_admm(bayer_noisy, conf)
%% call_demo_joint_demosaicing_denoising_admm.m
% TomHeaven, hanlin_tan@nudt.edu.cn, 2016.12.26
% Input: bayer_noisy - a (height * width) matrix denoting Bayer image with 
%        noise
%        conf - parameter structure for joint demosaicing denoising with
%        admm
% Output: RGB - a (height * width * 3) matrix for the result RGB image

addpath('utils');
addpath('admm');

conf.lambda_data = 1.1;   % data weight
conf.lambda_tv = 0.0065;  % tv weight
conf.lambda_cc = 0.06;    % cross-channel weight
conf.lambda_dm = 0.06;    % matlab demosaic weight
conf.beta = [0.25 0.25 0.11]; % channel weights
conf.rho = 0.5;           % soft thresholding weight

[height, ~, ~]= size(bayer_noisy);

I_raw_noisy = bayer_noisy(:);

% linear transform to 0 - 1
x = double(I_raw_noisy) / 255.0;

% get bayer image vector b
b = x(:);

% minimize using ADMM
x= joint_demosaicing_denoising_admm(conf.A, b, height, conf);

% get image X from vector x
X = vec2rgb(x, size(bayer_noisy, 2));

% linear tansform back to 0 - 255
RGB = uint8(X * 255);