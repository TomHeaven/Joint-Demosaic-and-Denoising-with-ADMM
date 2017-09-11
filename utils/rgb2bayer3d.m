function bayer3d = rgb2bayer3d(rgb)
%% transform a rgb image to bayer 3d image by setting pixels to zero 
% Input: rgb - height * width * 3, a rgb image
% Output: bayer3d - height * width * 3, a bayer image in 3d form

filter = zeros(size(rgb));
% R
filter(1:2:end, 1:2:end, 1) = 1;
% G
filter(2:2:end, 1:2:end, 2) = 1;
filter(1:2:end, 2:2:end, 2) = 1;
% B
filter(2:2:end, 2:2:end, 3) = 1;

bayer3d = double(rgb) .* filter;

end