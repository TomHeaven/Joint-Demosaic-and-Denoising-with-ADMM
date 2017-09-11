function x = joint_demosaicing_denoising_admm(A, b, height, conf)
%% Solve joint ISP minimization via ADMM of multile terms
%
%   minimize_x  lambda*||A*x - b||_2^2 + lambda_tv * sum_i (|x_{i+1} - x_i| +
%   | x{i+height} - x_i |) +  lambda_bm3d * bm3d(x) + lambda_cc * cc(x) +
%   lambda_dm * demosaic(x)
%
% where b in R^n is the Bayer image vector.
% The solution is returned in the vector x.
%
% Input: 
%        A - add mosaic matrix
%        b - a (height * width) vector denoting Bayer image with noise
%        conf - parameter structure for joint demosaicing denoising with
%        admm
% Output: x - a (height * width * 3) vector for the result RGB image

% By TomHeaven, hanlin_tan@nudt.edu.cn, 2016.11.09


%% config
MAX_ITER = 50;
EPS = 0.01;
nItem = 6;

lambda_data = conf.lambda_data;
lambda_tv = conf.lambda_tv;
lambda_cc = conf.lambda_cc;
lambda_dm = conf.lambda_dm;
beta = conf.beta;
rho = conf.rho;
sigma = conf.sigma;

%% init

addpath('admm/BM3D');
t_start = tic;

mu = rho;
if sigma < 1e-8 % sigma cannot be zero
    sigma = 0.1;
end

% difference matrix
n = 3 * length(b);
width = length(b) / height;

e = ones(n,1);
D1 = spdiags([e -e], [0 1], n,n);
D2 = spdiags([e -e], [0 height], n,n);
% init H
H = cell(nItem, 1);
H{1} = A;
H{2} = D1; 
H{3} = D2;
H{4} = spdiags(e, 0, n, n);

z0 = [b; b; b]; % need udpate
H{5} = crossChannelMatrix(z0, beta); 
H{6} = H{4};

b_gray = vec2gray(b, width); 
b_dm = demosaic(uint8(b_gray * 255), 'rggb');
% b_dm = call_onestep_ap(b_gray * 255);
b_dm = double(b_dm) ./ 255;
b_dm = rgb2vec(b_dm);


% init zeta, u, d, gamma
zeta = cell(nItem, 2); % 2 for k+1, 1 for k
u = cell(nItem, 2);
d = cell(nItem, 2);
gamma = cell(2, 1);
z = cell(2, 1);

vec0 = zeros(n, 1);
mat0 = sparse(n, n);
for k = 1 : 2
    for j = 1 : nItem
        zeta{j, k} = zeros(size(H{j}, 1), 1);   
        u{j,k} =  zeta{j, k};  
        d{j,k} =  zeta{j, k}; 
    end
    gamma{k} = vec0;
    %z{k} = b;
end


%% Iterations
for k = 1:MAX_ITER
    % update zeta
    for j = 1 : nItem
        zeta{j, 1} = u{j, 1} + d{j, 1};  
    end
    % update gamma
    gamma{1} = vec0;
    for j = 1 : nItem
        gamma{1} = gamma{1} + H{j}' * zeta{j, 1};
    end
    % update z
    sumH = sparse(mat0);
    for j = 1 : nItem
        sumH = sumH + H{j}'*H{j};
    end

    z{2} = sumH \ gamma{1};
  
    % update v_k and u_{k+1}
    
    % j = 1 data term
    nu = H{1} * z{2} - d{1, 1};
    u{1, 2} = (mu * nu + lambda_data*b) / (mu + lambda_data);
    % j = 2 vertical TV term
    nu = H{2} * z{2} - d{2, 1};
    u{2, 2} = shrinkage(nu, lambda_tv/rho);
    % j = 3 horizental TV term
    nu = H{3} * z{2} - d{3, 1};
    u{3, 2} = shrinkage(nu, lambda_tv/rho);
    % j = 4 CBM3D for denoising
    nu = H{4} * z{2} - d{4, 1};
    zRGB = vec2rgb(nu, width);
    [~, yRGB_est] = CBM3D(1, zRGB, sigma);
    u{4,2} = rgb2vec(yRGB_est);
    % j = 5 cross channel term
    H{5} = crossChannelMatrix(z{2}, beta);
    nu = H{5} * z{2} - d{5, 1};
    u{5, 2} = shrinkage(nu, lambda_cc/rho);
    % j = 6 demosaic
    nu = H{6} * z{2} - d{6, 1};
    u{6, 2} = (mu * nu + lambda_dm*b_dm) / (mu + lambda_dm);
    % update d_k
    err = 0;
    for j = 1 : nItem
        d{j, 2} = d{j, 1} - (H{j} * z{2} - u{j, 2}); 
        err = err + norm(d{j,2} - d{j, 1});
    end
    
    fprintf('Iter %d, err = %f\n', k, err);
    % check stopping criteria
    if err < EPS
        break;
    end
    
    % update k to k+1
    for j = 1 : nItem
        d{j,1} = d{j, 2};
        u{j,1} = u{j, 2};
        zeta{j, 1} = zeta{j, 2};
    end
    gamma{1} = gamma{2};
    z{1} = z{2}; 
end
%% return result
x = z{2};
% count time
toc(t_start);

end % end function

% soft threshold operator
function y = shrinkage(a, kappa)
    y = max(0, a-kappa) - max(0, -a-kappa);
end