function compareAlgorithms()
%% compareAlgorithms.m
% This file compares our method with FlexISP and Deep Joint.
%
% By TomHeaven, hanlin_tan@nudt.edu.cn, 2016.12.26

close all;
%% environment
addpath('utils');
addpath('wrappers');
addpath('admm');
addpath('deep_admm');
% caffe python interface path required by Deep Joint method
setenv('PYTHONPATH', '/Users/tomheaven/Documents/caffe-master/python:/home/tomheaven/文档/caffe/python'); 
% for ubuntu
setenv('LD_LIBRARY_PATH', ''); 

%% control
conf.debug = false;          % enable debug mode
conf.enableADMM = true;      % enable admm method
conf.enableFlexISP = true;   % enable flex isp method
conf.enableDeepJoint = true; % enable deep joint method
conf.cpuorgpu = 'gpu';    % 'gpu' or 'cpu', wheter to use gpu in deep joint
conf.tmpdir = '~/tmp/';   % a directory to store temporary output

%% config
conf.sigma = 0;       % noise levels: 0, 5, 15, 25
% test dataset
conf.dataset = 'Kodak';
%conf.dataset ='McM';

%% processing 
start = 1;
if strcmp(conf.dataset, 'McM')
    dataDir = 'data/McM/';  fileExt = '*.tif';
else
    dataDir = 'data/kodak/'; fileExt = '*.png'; 
end

files = dir(fullfile(dataDir, fileExt)); 

% sort filenames in numerical order
filenames = cell(length(files), 1); 
for i = 1 : length(files) 
    fileName = strcat(dataDir, files(i,1).name); 
    filenames{i} = fileName; 
end
filenames = natsortfiles(filenames); 

if conf.debug
    len = 1;
else
    len = length(filenames);
end

savepath = sprintf('res/%s/compareResults_%s_sigma%d.mat', conf.dataset, conf.dataset, conf.sigma);
for i=start:len
    fileName = filenames{i};
    fprintf('\nProcessing image %d, path = %s\n', i, fileName);
    I_groundtruth = imread(fileName);
    
    if conf.debug
       % I_groundtruth = imresize(I_groundtruth, 0.3); % 0.25 for kodak
       %I_groundtruth = imresize(I_groundtruth, [768 512]);
    end
    res = compare(I_groundtruth, conf, i);
    res.groundtruth_path = fileName;
    results{i} = res;
    save(savepath, 'results');
end

% analyze results
analyze(conf.dataset, conf.sigma, 3); 

end % func compareAlgorithms

function res = compare(I_groundtruth, conf, i)
%% compare algorithms using one input groudtruth
width = size(I_groundtruth, 2);
height = size(I_groundtruth, 1);
% get add-mosaic matrix A
    % caution: with frist, then height!!!
A = addMosaicMatrix(width, height);
conf.A = A;
% get image vector with mosaic
%x_raw_input = A * double(I_groundtruth(:));
% vector to rgb image
%I_raw_input = vec2gray(x_raw_input, width);

I_raw_input = sum(rgb2bayer3d(I_groundtruth), 3);
% add noise
rng(0);
bayer_noisy  =  uint8(I_raw_input + conf.sigma * randn(size(I_raw_input)) + 0.5);

% for algorithms that load image from file
noisypath = sprintf('data/noisy%d_sigma%d.tiff', i, conf.sigma);
imwrite(bayer_noisy, noisypath);

res.groundtruth = I_groundtruth;
res.bayer_noisy = bayer_noisy;
res.conf = conf;
res.noisypath = noisypath;


%bayer_noisy = im2double(bayer_noisy);


% 1. ADMM method
if conf.enableADMM
    s = tic;
    RGB = call_joint_demosaicing_denoising_admm(bayer_noisy, conf);
    time = toc(s);
    res.admm.psnr = comppsnr_rgb(I_groundtruth, RGB);
    res.admm.RGB = RGB;
    res.admm.time = time;
end

% 2. FlexISP
if conf.enableFlexISP
    s = tic;
    RGB = call_flexisp(bayer_noisy, I_groundtruth, conf);
    time = toc(s);
    res.flexisp.psnr = comppsnr_rgb(I_groundtruth, RGB);
    res.flexisp.RGB = RGB;
    res.flexisp.time = time;
end


% 3. Deep joint
if conf.enableDeepJoint
    s = tic;
    
    if conf.sigma > 19.1
        noiselevel = 19.1 / 255.0;
    else
        noiselevel = conf.sigma / 255.0;
    end
    commandStr = sprintf('python demosaicnet/demosaick --input %s --output data --model demosaicnet/pretrained_models/bayer_noise --noise %.2f --%s --offset_x %d', ...
        noisypath, noiselevel, conf.cpuorgpu, 1);
    [status, commandOut] = system(commandStr);
    time = toc(s);
    if status~=0
        fprintf('command failed! Log is:\n%s\n', commandOut);
        res.dj.failed = true;
        res.dj.log = commandOut;
    else
        RGB = imread('data/dj_output.jpg');
        res.dj.psnr = comppsnr_rgb(I_groundtruth, RGB);
        res.dj.RGB = RGB;
        res.dj.time = time;
        res.dj.failed = false;
        delete('data/dj_output.jpg');
    end
end

% 4. Matlab demosaic
s = tic;
RGB = demosaic(bayer_noisy, 'RGGB');
time = toc(s);
res.matlab.RGB = RGB;
res.matlab.psnr = comppsnr_rgb(I_groundtruth, RGB);
res.matlab.time = time;

% delete temporary file
delete(noisypath);

end % func compare
