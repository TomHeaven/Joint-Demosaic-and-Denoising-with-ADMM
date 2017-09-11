function X = vec2rgb(x, width)
% Get rgb image matrix from vector x = [r, g, b]'
%% Input:
%    x - (3MN)*1 vector
%    width - the width of the original rgb image
%  Output:
%    X - M*N*3 image matrix
% By TomHeaven, hanlin_tan@nudt.edu.cn @ 2016.10.17

n = int32(length(x)/3);
N = int32(width);
M = int32(n / width);
% get r,g,b
r = x(1:n);
g = x(n+1:2*n);
b = x(2*n+1 : 3*n);
% reshape
R = reshape(r, [M N]);
G = reshape(g, [M N]);
B = reshape(b, [M N]);
% get X
X = zeros(M,N,3);
X(:,:,1) = R;
X(:,:,2) = G;
X(:,:,3) = B;

end