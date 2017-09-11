function X = vec2rgb(x, width)
% Get rgb image matrix from vector x = [r, g, b]'
%% Input:
%    x - (MN)*1 vector
%    width - the width of the original rgb image
%  Output:
%    X - M*N image matrix
% By TomHeaven, hanlin_tan@nudt.edu.cn @ 2016.10.17

n = int32(length(x));
N = int32(width);
M = int32(n / width);

% reshape
X = reshape(x, [M N]);


end