addpath('../utils');

dims = [512 768];

%I = reshape(1:dims(1)*dims(2)*3, [dims(1), dims(2), 3]);

I = imread('data/kodak/kodim01.png');
I = double(I);

i = I(:);

A = addMosaicMatrix(size(I, 2), size(I, 1));
b = A * i;

B = vec2gray(b, size(I, 2));
B2 = sum(rgb2bayer3d(I), 3);

figure; imshow(uint8(I)); title('I')
figure; imshow(uint8(B2)); title('B2')

disp(norm(B - B2))