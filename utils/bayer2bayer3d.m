function bayer3d = bayer2bayer3d(bayer)
%% transform a rgb image to bayer 3d image by setting pixels to zero 
% Input: bayer - height * width, a rgb image
% Output: bayer3d - height * width * 3, a bayer image in 3d form

RGB = zeros(size(bayer, 1), size(bayer, 2), 3);
for i = 1:3
    RGB(:,:,i) = bayer;
end


bayer3d = rgb2bayer3d(RGB);

end