function [CR, CG, CB] = buildCRGB(width, height)
% Get rgb image matrix from vector x = [r, g, b]'
%% Input:
%    width - the width of the original rgb image
%    height - the height of the original rgb image
%  Output:
%    CR, CG, CB - height*width lable image matrix to filter R,G,B channels 
%    from Bayer image
% By TomHeaven, hanlin_tan@nudt.edu.cn @ 2016.10.17

% build CR,CG,CB matrices 
CR = toeplitz( mod( 1:height, 2), mod( 1:width, 2) ); 
CR(2:2:end, :) = 0;
CG = toeplitz( mod( 0:height - 1, 2), mod( 0:width - 1, 2) );
CB = toeplitz( mod( 1:height, 2), mod( 1:width, 2) ); 
CB(1:2:end - 1, :) = 0;

end