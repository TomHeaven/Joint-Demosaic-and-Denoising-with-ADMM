function x = rgb2vec(X)
% Transform rgb image matrix to vector x = [r, g, b]'
%% Input:
%    X - M*N*3 image matrix
%  Output:
%    x - (3MN)*1 vector
% By TomHeaven, hanlin_tan@nudt.edu.cn @ 2016.10.17

R = X(:,:,1);
G = X(:,:,2);
B = X(:,:,3);

% transform into column prior vectors
r = R(:); 
g = G(:);
b = B(:);

x = [ r; g; b];


end