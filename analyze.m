function analyze(dataset, sigma, nAlgos)
%% analyze and visualize results of compareAlgorithms
% By TomHeaven, hanlin_tan@nudt.edu.cn, 2016.12.27

%close all;

%% settings
fontSize = 20; % figure font size
enableExportResultImages = true;
% combine results with previous results stored in 'compareResults_nonAtomic.mat'
enableUpdateResultMat = false;
% analyze this number of algorithms
% nAlgos = 3;

if enableUpdateResultMat
    load('compareResults_old.mat');
    results0 = results;
end

savepath = sprintf('res/%s/compareResults_%s_sigma%d.mat', dataset, dataset, sigma);
load(savepath);
% sort results according to filenames
results = sortResultsFileNames(results);

len = length(results);

if ~exist('nAlgos', 'var')
    nAlgos = 5;
end

if enableUpdateResultMat
    %len = 2; % debug

    for i = 1 : len
        if ~isfield(results{i}, 'dadmm') && isfield(results0{i}, 'dadmm') == 1
            results{i}.dadmm = results0{i}.dadmm;
        end
        if ~isfield(results{i}, 'admm') && isfield(results0{i}, 'admm') == 1
            results{i}.admm = results0{i}.admm;
        end
        if ~isfield(results{i}, 'flexisp') && isfield(results0{i}, 'flexisp') == 1
            results{i}.flexisp = results0{i}.flexisp;
        end
        if ~isfield(results{i}, 'dj') && isfield(results0{i}, 'dj') == 1
            results{i}.dj = results0{i}.dj;
        end
        if ~isfield(results{i}, 'fcn') && isfield(results0{i}, 'fcn') == 1
            results{i}.fcn = results0{i}.fcn;
        end
        if ~isfield(results{i}, 'fcne') && isfield(results0{i}, 'fcne') == 1
            results{i}.fcne = results0{i}.fcne;
        end
        results{i}.groundtruth = results0{i}.groundtruth;
        results{i}.bayer_noisy = results0{i}.bayer_noisy;
        results{i}.conf = results0{i}.conf;
    end
    save(savepath, 'results');
end

psnr = zeros(len, nAlgos);
time = zeros(len, nAlgos);

if results{1}.conf.debug
    len = 1;
end

for i = 1 : len
    if isfield(results{i}, 'admm') == 1
        psnr(i, 1) = results{i}.admm.psnr;
        time(i, 1) = results{i}.admm.time;
    end
    if isfield(results{i}, 'flexisp') == 1
        psnr(i, 2) = results{i}.flexisp.psnr;
        time(i, 2) = results{i}.flexisp.time;
    end
    if isfield(results{i}, 'dj') == 1
        psnr(i, 3) = results{i}.dj.psnr;
        time(i, 3) = results{i}.dj.time;
    end
    if nAlgos > 3
        if isfield(results{i}, 'fcn') == 1
            psnr(i, 4) = results{i}.fcn.psnr;
            time(i, 4) = results{i}.fcn.time;
        end
        if isfield(results{i}, 'fcne') == 1
            psnr(i, 5) = results{i}.fcne.psnr;
            time(i, 5) = results{i}.fcne.time;
        end
        if isfield(results{i}, 'dadmm') == 1
        psnr(i, 6) = results{i}.dadmm.psnr;
        time(i, 6) = results{i}.dadmm.time;
    end
    end
end

ave_pnsr =  sum(psnr, 1) / len;

fprintf('avg psnr for admm = %f\n',ave_pnsr(1) );
fprintf('avg psnr for flexisp = %f\n', ave_pnsr(2));
fprintf('avg psnr for dj = %f\n', ave_pnsr(3));

if nAlgos > 3
    if isfield(results{i}, 'fcne') == 1
        fprintf('avg psnr for fcn = %f\n', ave_pnsr(4));
    end
    if isfield(results{i}, 'fcne') == 1
        fprintf('avg psnr for fcn_ensemble = %f\n', ave_pnsr(5));
    end
    fprintf('avg psnr for deep_admm = %f\n', ave_pnsr(6));
end



h1 = figure;
bar(psnr);
if nAlgos == 5
    legend('ADMM', 'FlexISP', 'DeepJoint', 'FCN', 'FCN-Ensemble', 'DeepADMM', 'Location', 'southeast');
else
    legend('ADMM', 'FlexISP', 'DeepJoint', 'Location', 'southeast');
end
title('PSNR comparison');
xlabel('Image No.');
ylabel('dB');
set(gca, 'FontSize', fontSize);
saveas(h1, sprintf('res/%s/psnr_comparison_%s_sigma%d_algos_%d.jpg', dataset, dataset, sigma, nAlgos));


h2 = figure;
bar(time);
if nAlgos == 5
    legend('ADMM', 'FlexISP', 'DeepJoint', 'FCN', 'FCN-Ensemble', 'DeepADMM', 'Location', 'east');
else
    legend('ADMM', 'FlexISP', 'DeepJoint', 'Location', 'east');
end

title('Running time comparison');
xlabel('Image No.');
ylabel('sec');
set(gca, 'FontSize', fontSize);
saveas(h2, sprintf('res/%s/time_comparison_%s_sigma%d_algos_%d.jpg', dataset, dataset, sigma, nAlgos));

if enableExportResultImages
    mkdir(sprintf('res/%s/sigma%d', dataset, sigma));
    for i = 1 : len
        if isfield(results{i}, 'admm') == 1
            imwrite(results{i}.admm.RGB, sprintf('res/%s/sigma%d/im%d_admm.bmp', dataset, sigma, i));
        end
        if isfield(results{i}, 'flexisp') == 1
            imwrite(results{i}.flexisp.RGB, sprintf('res/%s/sigma%d/im%d_flexisp.bmp', dataset, sigma, i));
        end
        if isfield(results{i}, 'dj') == 1
            imwrite(results{i}.dj.RGB, sprintf('res/%s/sigma%d/im%d_deepjoint.bmp', dataset, sigma, i));
        end
        if isfield(results{i}, 'fcn') == 1
            imwrite(results{i}.fcn.RGB, sprintf('res/%s/sigma%d/im%d_fcn.bmp', dataset, sigma, i));
        end
        if isfield(results{i}, 'fcne') == 1
            imwrite(results{i}.fcne.RGB, sprintf('res/%s/sigma%d/im%d_fcne.bmp', dataset, sigma, i));
        end
        if isfield(results{i}, 'dadmm') == 1
            imwrite(results{i}.dadmm.RGB, sprintf('res/%s/sigma%d/im%d_dadmm.bmp', dataset, sigma, i));
        end
        imwrite(uint8(results{i}.bayer_noisy), sprintf('res/%s/sigma%d/im%d_bayer_noisy.bmp', dataset, sigma, i));
        imwrite(results{i}.groundtruth, sprintf('res/%s/sigma%d/im%d_groundtruth.bmp', dataset, sigma, i));
        imwrite(results{i}.matlab.RGB, sprintf('res/%s/sigma%d/im%d_matlab.bmp', dataset, sigma, i));
    end
end

for i = 1 : len
    res = results{i};
    sz = size(res.groundtruth);
    %    name = sscanf(results{i}.noisypath, '../data/%s');
    %   l = length(name);
    %   name = name(1:l - 5);
    
    % fprintf('%d & %s & $%d \\times %d $ & %.2f / %.3f & %.2f / %.2f & %.2f / %.2f & %.2f / %.2f & %.2f / %.2f \\\\ \n', ...
    %     i, name, sz(1), sz(2), res.admm.psnr, res.admm.time, res.flexisp.psnr, res.flexisp.time, ...
    %     res.dj.psnr, res.dj.time, res.fcn.psnr, res.fcn.time, res.fcne.psnr, res.fcne.time);
end

save('analysis.mat', 'psnr', 'time', 'ave_pnsr');
end


function results = sortResultsFileNames(results)
    addpath('utils');
    
    len = length(results);
    paths = cell(len, 1);
    for i = 1 : len
        paths{i} = results{i}.groundtruth_path;
    end
    paths = natsortfiles(paths);
    for i = 1 : len
        for j = 1 : len
            if strcmp(results{j}.groundtruth_path, paths{i}) == 1 && i~= j
                t = results{i};
                results{i} = results{j};
                results{j} = t;
                break;
            end
        end
    end
end