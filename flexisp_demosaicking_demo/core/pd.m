% Copyright 2014 NVIDIA Corporation.  All Rights Reserved
% This software is for provided for non-commercial and research purposes only
% Paper: FlexISP: A Flexible Camera Image Processing Framework, SIGGRAPH
% Asia 2014
% WebSite: https://research.nvidia.com/publication/flexisp-flexible-camera-image-processing-framework
%%%%%%%%%%%%%%%%%%
% Special Thanks: 
% Gabriel Peyre. He has made great contribution to the scientific community on numerical analysis of image processing
% Some of the code is derived from here: https://www.ceremade.dauphine.fr/~peyre/numerical-tour/tours/

function db_chs = pd(channels, channels_0, ...
                     w_channels, w_regularization, w_sigma_scale, max_it,    ...
                     tol, output_folder, verbose, I_sharp)

%Check for channel sanity
if ~isempty(channels_0) && (length(channels_0) ~= length(channels))
  error('Initial channels do not match channels.\n')
end

%Debug options
local_iterate_fig = [];
if strncmpi(verbose,'all',3)
    verbose = 'all';
    %Figures
    local_iterate_fig = figure();
elseif strncmpi(verbose,'brief',5)
    verbose = 'brief';
else
    verbose = 'none';
end           

%Initialize all channels
if isempty(channels_0)
    channels_0 = channels;
end

db_chs = channels; %Empty initialization
for ch = 1:length(channels)
    
    %Set to the initial iterate if neccessary
    if ~isempty(db_chs(ch).K)
        db_chs(ch).Image = channels_0(ch).Image;
    elseif isequal(size(db_chs(ch).K, 1), size(db_chs(ch).K, 2)) || mod(size(db_chs(ch).K, 1),2) == 0
        error('Blur size is not square or odd'); %Blur size check
    end
                    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Do startup minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Do startup minimization
if strcmp(verbose, 'brief') || strcmp(verbose, 'all')
    fprintf('\n### Startup iteration ###\n')
end

%Get current channel and weights
w_res_curr =   w_channels(:, 1);  
w_cross_curr = w_channels(:, 2:end);

% edgetaper to better handle circular boundary conditions
%ks = size(channels(ch_opt).K,1);
for ch = 1:length(channels)

    channels(ch).Image = channels(ch).Image + 1.0;
    db_chs(ch).Image = db_chs(ch).Image + 1.0;
end

%Do residual pd deconvolution
[db_chs] = pd_channel_deconv(channels, db_chs, [], ...
                      w_res_curr, w_cross_curr, w_regularization, w_sigma_scale, ...
                      max_it, tol, verbose, output_folder, local_iterate_fig, I_sharp);  

for ch = 1:length(channels)

    channels(ch).Image = channels(ch).Image - 1.0;
    db_chs(ch).Image = db_chs(ch).Image - 1.0;
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Main routine. Please refer to the supplemental material Algorithm 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ db_chs ] = pd_channel_deconv(channels, db_chs, x_0, ...
                    lambda_residual_cap, lambda_cross_ch, lambda_prior, sigma_scale, ...
                    max_it, tol, verbose, output_folder, iterate_fig, I_sharp)
                
    %Generate output folder
    if ~isempty( output_folder )
        mkdir( output_folder, 'iterations' );
    end
    
    %Convert to an image for simplicity
    chImg =   cat( 3, channels(1).Image, channels(2).Image, channels(3).Image );
    dbchImg = cat( 3, db_chs(1).Image,     db_chs(2).Image,   db_chs(3).Image );
    
    %Shortcut for the TV norm.
    Amplitude = @(u)sqrt(u.^2);
    
    %Generate both captures for the HDR capture
    MtB_cap = observationMat( chImg, channels, -1  );
    
    %Lambda TV
    lambda_tv = 0.0; %5;
    
    %Prox operator    
    %L1 norm
    ProxL1S = @(u,gamma) u ./ max(1, Amplitude(u));

    %Prox BM3D
    bm3d_fig = [];
    if strcmp(verbose, 'all')
        bm3d_fig = figure();
    end
    ProxBM3DS = @(v,gamma, it) v - gamma * BM3D_proc(sqrt(1/gamma), v / gamma, lambda_prior, sigma_scale, bm3d_fig, output_folder, it);
    
    %ProxFS
    ProxFS = @(v, gamma, it) prox_wrap(v, gamma, it, ...
                                       lambda_prior, lambda_tv, lambda_cross_ch, ProxBM3DS, ProxL1S);

    %Fast data solve
    lambda_residual_cap = cat( 3, repmat( lambda_residual_cap(1), size(chImg,1), size(chImg,2) ), ...
                                  repmat( lambda_residual_cap(2), size(chImg,1), size(chImg,2) ), ...
                                  repmat( lambda_residual_cap(3), size(chImg,1), size(chImg,2) ) );
                              
    ProxG = @(f,tau) solve_cg(f, channels, MtB_cap, tau, lambda_residual_cap );
    
    %Define A and AS
    A = @(v) Kmult_channels(v, db_chs, lambda_cross_ch, lambda_prior, lambda_tv);
    AS = @(v) KSmult_channels(v, db_chs, lambda_cross_ch, lambda_prior, lambda_tv );
    
    %Minimization weights
    L = compute_operator_norm(A, AS, size(chImg));
    sigma = 27.0;
    
    if sigma_scale < 0
        sigma = 10;
    end
    
    tau = 0.9/(sigma*L^2);
    theta = 1.0;
    
    %Set initial iterate
    if isempty(x_0)
        x_0 = dbchImg;
    end
    
    f = x_0;
    g = A(f);
    f1 = f;
        
    %Primal-Dual Iterative solver
    for i = 1:max_it
        
        fold = f;
        
        %Get prox input
        gp = A(f1);
        for ch = 1:length(gp)
            gp{ch} = g{ch} + sigma * gp{ch};
        end
        
        g = ProxFS( gp, sigma, i);
        f = ProxG( f - tau * AS(g), tau);
        f1 = f + theta * (f-fold);
        
        %Compute PSNR
        PSNR = 0;
        if ~isempty(I_sharp)
            
            %Compute PSNR:
            f_test = max(min(f1 - 1,1),0);
            
            psnr_pad = round(size(channels(1).K,1) * 1.0);
            I_diff = I_sharp(psnr_pad + 1:end - psnr_pad, psnr_pad + 1:end - psnr_pad,:) - f_test(psnr_pad + 1:end - psnr_pad, psnr_pad + 1:end - psnr_pad,:);
            MSE = 1/size(I_diff(:),1)*(norm(I_diff(:), 2)^2);
            if MSE > eps
                PSNR = 10*log10(1/MSE);
            else
                PSNR = Inf;
            end
            
        end
        
        %Save
        if ~isempty( output_folder )
            imwrite(f1 - 1, sprintf( '%s/iterations/iter%03d_PSNR%2.2f.png', output_folder, i, PSNR ), 'Bitdepth', 16 );
        end
        
        diff = f - fold ;
        f_comp = f;
        if strcmp(verbose, 'brief') || strcmp( verbose, 'all')
            fprintf('Iter %d, diff %5.5g, PSNR %2.2f\n', i, norm(diff(:),2)/ norm(f_comp(:),2), PSNR)
        end
        if norm(diff(:),2)/ norm(f_comp(:),2) < tol
            break;
        end
    end
    
    %Threshold if necessary
    f1( f1 < 1) = 1;
    
    %Copy back
    db_chs(1).Image = f1(:,:,1);
    db_chs(2).Image = f1(:,:,2);
    db_chs(3).Image = f1(:,:,3);
    
return;

function x = solve_cg(f, channels, MtB_cap, tau, lambda_residual_cap )

    %Solves Ax = b with
    % A = (lambda_residual_cap*tau* M'* M + eye ) and b = tau * lambda_residual_cap * M' * b + f
    % Please refer to eq(7) in the paper
    
    %Right matrix side
    b = tau * 2*lambda_residual_cap .* MtB_cap + f;
    
    %Matrix
    Mfun = @(x) reshape( tau * 2*lambda_residual_cap .* observationMat( reshape(x,size(f)) , channels, 0  ) + reshape(x,size(f)), [], 1 );  

    %Solve
    warning off;
    [x, ~, ~, ~] = pcg( Mfun, b(:), 1e-12, 100, [], [], f(:) );
    warning on;
    
    %Reshape
    x = reshape( x, size(f) );

return;


function prox = prox_wrap(v, gamma, it, ...
                          lambda_prior, lambda_tv, lambda_cross_ch, ProxBM3DS, ProxL1S)
    
    %Resulting prox
    prox = v;
    
    i = 0;
    if lambda_prior > eps
        
        %Inc
        i = i + 1;
        
        %Get image
        v_img = [];
        for ch = 1:length(v)
            v_img = cat( 3, v_img, v{ch}(:,:, i) );
        end

        %Get prox
        prox_bm3d = ProxBM3DS( v_img, gamma, it);
        
        %Insert
        for ch = 1:length(v)
            prox{ch}(:,:, i) = prox_bm3d(:,:, ch);
        end

    end

    %Compute tv terms
    if lambda_tv > eps

        %Gather
        i = i + 1;
        
        %Apply
        for ch = 1:length(v)
            prox{ch}(:,:, i:i + 1) = ProxL1S( v{ch}(:,:, i:i + 1), gamma);
        end
        
        %Add i
        i = i + 1;

    end

    %Cross channel prior 
    %Gather
    i = i + 1;

    %Apply
    for ch = 1:length(v)

        %Check for applicable end
        if sum( lambda_cross_ch(ch,:) ) <= eps
            continue;
        end

        prox{ch}(:,:,i:end) = ProxL1S( v{ch}(:,:,i:end), gamma);
    end
    
    %Default fill as identity
    if isempty(prox)
        warning('Prox seems to be empty, using identity.')
        prox = v;
    end
    
return;

function Kmultf = Kmult_channels(f, db_chs, lambda_cross_ch, lambda_prior, lambda_tv )
    
    %Gather
    for ch = 1:size(f,3)
        Kmultf{ch} = Kmult( f(:,:,ch), ch, db_chs, lambda_cross_ch(ch,:), lambda_prior, lambda_tv);
    end

return;

function Kmultf = Kmult(f, ch, db_chs, lambda_cross_ch, lambda_prior, lambda_tv )

%Result
Kmultf = [];

% derivative filters
dxf=[-1 1];
dyf=[-1;1];

%Compute tv terms
if lambda_prior > eps

    %Gather
    Kmultf = lambda_prior * f;
    
end

if lambda_tv > eps
    
    %Operators for anisotropic TV
    S_tv = lambda_tv * grad(f);
    Kmultf = cat(3,Kmultf, S_tv);
    
end

%Cross-Terms for all adjacent channels
for adj_ch = 1:length(db_chs)

    %Continue for current channel and zero channels
    if adj_ch == ch || lambda_cross_ch(adj_ch) < eps()
        continue;
    end
    adjChImg = db_chs(adj_ch).Image; %Curr cross channel
    
    %Compute cross terms
    diag_term = imconv(adjChImg, fliplr(flipud(dxf)), 'full');
    diag_term = diag_term(:, 2:end) .* f;
    conv_term = imconv(f,fliplr(flipud(dxf)), 'full');
    Sxf = (lambda_cross_ch(adj_ch) * 0.5) * ( adjChImg .* conv_term(:, 2:end) - diag_term );
    
    diag_term = imconv(adjChImg, fliplr(flipud(dyf)), 'full');
    diag_term = diag_term(2:end, :) .* f;
    conv_term = imconv(f,fliplr(flipud(dyf)), 'full');
    Syf = (lambda_cross_ch(adj_ch) * 0.5) * ( adjChImg .* conv_term(2:end, :) - diag_term );
    
    %Gather
    Kmultf = cat(3,Kmultf, Sxf, Syf);
end

%Default fill as identity
if isempty(Kmultf)
    Kmultf = f;
end

return;

function KSmultf = KSmult_channels(f, db_chs, lambda_cross_ch, lambda_prior, lambda_tv )
    
    %Resulting prox
    KSmultf = [];

    %Gather
    for ch = 1:length(f)
        KSmultf = cat(3, KSmultf, KSmult( f{ch}, ch, db_chs, lambda_cross_ch(ch,:), lambda_prior, lambda_tv) );
    end

return;

function KSmultf = KSmult(f, ch, db_chs, lambda_cross_ch, lambda_prior, lambda_tv )

%Result
KSmultf = zeros([size(f,1), size(f,2)]);

% derivative filters
dxf=[-1 1];
dyf=[-1;1];

%Compute tv terms
i = 0; %Image access term
if lambda_prior > eps

    i = i + 1;
    fpr = lambda_prior * f(:,:,i);

    % gather result
    KSmultf = KSmultf + fpr;
end

if lambda_tv > eps

    i = i + 1;
    f_tv = - lambda_tv * div(f(:,:,i:i+1));
    i = i + 1;
    
    % gather result
    KSmultf = KSmultf + f_tv;
end

%Cross-Terms for all adjacent channels
for adj_ch = 1:length(db_chs)
    
    %Continue for current channel and zero channels
    if adj_ch == ch || lambda_cross_ch(adj_ch) < eps()
        continue;
    end
    adjChImg = db_chs(adj_ch).Image; %Curr cross channel
    
    %Compute cross terms
    i = i + 1;
    f(:,:,i) = (lambda_cross_ch(adj_ch) * 0.5) * f(:,:,i);
    diag_term = imconv(adjChImg, fliplr(flipud(dxf)), 'full');
    diag_term = diag_term(:, 2:end) .* f(:,:,i);
    conv_term = imconv(adjChImg .* f(:,:,i), dxf, 'full'); 
    Sxtf = ( conv_term(:,1:(end-1)) - diag_term );
    
    i = i + 1;
    f(:,:,i) = (lambda_cross_ch(adj_ch) * 0.5) * f(:,:,i);
    diag_term = imconv(adjChImg, fliplr(flipud(dyf)), 'full');
    diag_term = diag_term(2:end, :) .* f(:,:,i);
    conv_term = imconv(adjChImg .* f(:,:,i), dyf, 'full'); 
    Sytf = ( conv_term(1:(end-1),:) - diag_term ); 

    %Gather result
    KSmultf = KSmultf + Sxtf + Sytf;
end

%Default fill as identity
if (lambda_prior <= eps) && (lambda_tv <= eps) && (sum(lambda_cross_ch) <= eps)
    KSmultf = ones([size(f,1), size(f,2)]);
end

return;


