% Copyright 2014 NVIDIA Corporation.  All Rights Reserved
% This software is for provided for non-commercial and research purposes only
% Paper: FlexISP: A Flexible Camera Image Processing Framework, SIGGRAPH
% Asia 2014
% WebSite: https://research.nvidia.com/publication/flexisp-flexible-camera-image-processing-framework

function L = compute_operator_norm(A, AS, sx)

    % computes the operator norm for a linear operator AS on images with size sx, 
    % which is the square root of the largest eigenvector of AS*A.
    % http://mathworld.wolfram.com/OperatorNorm.html
    vec_size = prod(sx);

    %Compute largest eigenvalue (in this case arnoldi, since matlab
    %implementation faster than power iteration)
    opts.tol = 1.e-3;
    opts.maxit = 50;
    lambda_largest = eigs(@(x)ASAfun(x, A, AS, sx ), vec_size, 1,'LM', opts);
    L = sqrt(lambda_largest);

return;

function ASAx = ASAfun(x, A, AS, sx)
    x_img = reshape(x,sx);
    ASAx = AS(A(x_img));
    ASAx = ASAx(:);
return;