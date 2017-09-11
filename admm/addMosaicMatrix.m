% By TomHeaven, hanlin_tan@nudt.edu.cn @ 2016.10.17
function A = addMosaicMatrix(width, height)
%% Input:
 %   x - 3n * 1 vector of a color image, where n = width * height
 %   height - an integer, the original height of the image
 %   width - an integer, the original width of the image
%% Output:
 %   y -  n * 1 vector of a bayer image with RGGB format
 %   A - the add-mosaic matrix of RGGB format, a sparse matrix
 
 if mod(width, 2) == 1 || mod(height, 2) == 1
     error('Input image dimension must not be odd');
 end
 
n = width * height;
% build CR,CG,CB matrices 
[CR, CG, CB] = buildCRGB(width, height);

CR_hat = spdiags(CR(:), 0, n, n); % 列优先向量作为对角矩阵的元素
CG_hat = spdiags(CG(:), 0, n, n); 
CB_hat = spdiags(CB(:), 0, n, n);

A = [CR_hat CG_hat CB_hat];

end