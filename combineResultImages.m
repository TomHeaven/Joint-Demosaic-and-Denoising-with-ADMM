function combineResultImages(dataset)
%% combine result images of different noise levels into one-picture-multi-level images
% By TomHeaven, hanlin_tan@nudt.edu.cn, 2017.01.21

close all;

%% settings
% dataset name
%dataset = 'McMaster'; % 'Kodak' or 'McMaster'
% input image directory
dirpath = sprintf('res/%s', dataset);
% noise levels in the order of counter-clock
noiselevels = [15 25 0 5];
% input file extensions
fileExt = '*.bmp';
% epsilon
EPS = 0.0001;

%% work

lenNoise =  length(noiselevels);
% load image results
results = cell(lenNoise, 1);
for i = 1 : lenNoise
    subdirname = sprintf('sigma%d', noiselevels(i));
    dataDir = sprintf('%s/%s', dirpath, subdirname);
    files = dir(fullfile(dataDir, fileExt));
    len = size(files,1);
    images = cell(len, 1);
    if i == 1
        names = cell(len, 1);
    end
    for j = 1 : len
        filepath = sprintf('%s/%s', dataDir, files(j, 1).name);
        images{j} =imread(filepath);
        names{j} = files(j, 1).name;
    end
    results{i}.images = images;
    if i == 1
        results{1}.names = names;
    end
end

% combine results of different noise levels into one result
outputDir = sprintf('%s/combined', dirpath);
if exist(outputDir, 'file') ~= 7
    mkdir(outputDir);
end

resImages = cell(len, 1);
%len = 1;
for j = 1 : len
    [rows, cols, ~] = size(results{1}.images{j});
    imageMaps = createImageMaps(rows, cols, lenNoise, EPS);
    nChannels = size(results{1}.images{j}, 3);
    resImages{j} = uint8(zeros(rows, cols, nChannels));
    for i = 1 : lenNoise
        srcImg = results{i}.images{j};
        for k = 1 : nChannels
            resImages{j}(:,:,k) = resImages{j}(:,:,k) + imageMaps{i} .* srcImg(:,:,k);
        end
    end
end

% write results
for j = 1 : len
    imwrite(resImages{j}, sprintf('%s/%s_%s.jpg', outputDir, results{1}.names{j}, dataset));
end

end % function

function imageMaps = createImageMaps(rows, cols, n, EPS)
% create image maps to crop and combine images

imageMaps = cell(n, 1);
center.x = rows / 2;
center.y = cols / 2;

for i = 1 : n
    imageMap = zeros(rows, cols);
    startAngle = (i - 1) * (2*pi/n);
    endAngle = i * (2*pi/n);
    
    for r = 1 : rows
        for c = 1 : cols
            % move origin to image center
            x = r - center.x;
            y = c - center.y;
            angle = abs(atan(y / x));
            % find the correct quadrant
            if x < 0 && y >= 0
                angle = angle + pi/2;
            elseif x < 0 && y < 0
                    angle = angle + pi;
            elseif x > 0 && y < 0
                angle = angle + pi*3/2;
            end        
            % build map according to angle range
            if angle > startAngle + EPS && angle < endAngle - EPS
                imageMap(r, c) = 1;
            end
        end
    end
    imageMaps{i} = uint8(imageMap);
end

end % function createImageMaps