% test cross channel prior
% By TomHeaven, hanlin_tan@nudt.edu.cn, 2016.11.09
function CC = crossChannelMatrix(x, beta)
% Input:
%   x for [r; g; b] - 3n*1 vectors for channel r,g,b
%   beta - 3*1 vector for coefficents for rg, gb, br
% Output:
%   CC - cross channel prior matrix

len = length(x);
n = len / 3;

r = x(1:n);
g = x(n+1:2*n);
b = x(2*n+1: 3*n);

e = ones(n, 1);
% difference matrices
H1 = spdiags([e -e], [0 1], n, n);
H2 = spdiags([e -e], [0 n], n, n);
% diagnal matrices for r,g,b
D_r = spdiags(r, 0, n, n);
D_g = spdiags(g, 0, n, n);
D_b = spdiags(b, 0, n, n);
% temporal diagnal marices 
D_H1g = spdiags(H1 * g, 0, n, n);
D_H2g = spdiags(H2 * g, 0, n, n);
D_H1b = spdiags(H1 * b, 0, n, n);
D_H2b = spdiags(H2 * b, 0, n, n);
D_H1r = spdiags(H1 * r, 0, n, n);
D_H2r = spdiags(H2 * r, 0, n, n);
% zero matrix of 4n * n
Z = sparse(4*n, n);
% assign coefficients
beta_rg = beta(1);
beta_gb = beta(2);
beta_br = beta(3);

% build matrix R of 4n*n
R = [
  beta_rg * ( D_g * H1 - D_H1g );  
  beta_rg * ( D_g * H2 - D_H2g );
  beta_br * ( D_b * H1 - D_H1b );
  beta_br * ( D_b * H2 - D_H2b );
];
% build matrix G of 4n*n
G = [
  beta_rg * ( D_r * H1 - D_H1r );
  beta_rg * ( D_r * H2 - D_H2r );
  beta_gb * ( D_b * H1 - D_H1b );
  beta_gb * ( D_b * H1 - D_H2b );
];
% build matrix B of 4n*n
B = [
  beta_gb * ( D_g * H1 - D_H1g);
  beta_gb * ( D_g * H2 - D_H2g);
  beta_br * ( D_r * H1 - D_H1r);
  beta_br * ( D_r * H2 - D_H2r);
];
% build cross channel prior matrix
CC = [ R Z Z;
       Z G Z;
       Z Z B];

end
