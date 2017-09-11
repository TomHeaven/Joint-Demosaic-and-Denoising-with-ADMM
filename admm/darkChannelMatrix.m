% test dark channel prior
% By TomHeaven, hanlin_tan@nudt.edu.cn, 2016.11.15
function DC = darkChannelMatrix(x)
% Input:
%   x for [r; g; b] - 3n*1 vectors for channel r,g,b
% Output:
%   DC - dark channel prior matrix

len = length(x);
n = len / 3;

r = x(1:n);
g = x(n+1:2*n);
b = x(2*n+1: 3*n);

% temporal diagnal marices 
D_gb = spdiags(g .* b, 0, n, n);
D_rb = spdiags(r .* b, 0, n, n);
D_rg = spdiags(r .* g, 0, n, n);

% zero matrix of 4n * n
Z = sparse(n, n);

% build dark channel prior matrix
DC = [ D_gb Z     Z;
       Z    D_rb  Z;
       Z    Z     D_rg ];

end
