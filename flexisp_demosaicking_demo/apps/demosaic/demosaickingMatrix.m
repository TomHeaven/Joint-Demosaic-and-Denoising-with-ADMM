% Copyright 2014 NVIDIA Corporation.  All Rights Reserved
% This software is for provided for non-commercial and research purposes only
% Paper: FlexISP: A Flexible Camera Image Processing Framework, SIGGRAPH
% Asia 2014
% WebSite: https://research.nvidia.com/publication/flexisp-flexible-camera-image-processing-framework

function [ result ] = demosaickingMatrix( x, ch, sizeX, flag  )

    %Joint HDR capture subsampling and demosaicking matrix
    
    %Channel mapping
    % r -> 1
    % g -> 2
    % b -> 3

    % check for x divisible by two
    if flag >= 0 && ( mod( size(x,1), 2 ) ~= 0 || mod( size(x,2), 2 ) ~= 0 )
        error('The image to mosaick is not divisible by two, which it should.')
    end
    

    %Build G matrix
    G_Mos = toeplitz( mod( 0:sizeX(1) - 1, 2), mod( 0:sizeX(2) - 1, 2) );

    R_Mos = toeplitz( mod( 1:sizeX(1), 2), mod( 1:sizeX(2), 2) ); 
    R_Mos(2:2:end, :) = 0;

    B_Mos = toeplitz( mod( 1:sizeX(1), 2), mod( 1:sizeX(2), 2) ); 
    B_Mos(1:2:end - 1, :) = 0;

    %Logical
    G_Mos = logical(G_Mos);
    R_Mos = logical(R_Mos);
    B_Mos = logical(B_Mos);
    
    % Diagonal matrix
    if flag > 0 || flag < 0 || flag == 0      
        
        % A'A
        result = zeros( sizeX );
        if ch == 1 
            %R-channel
            result(R_Mos) = reshape(x(R_Mos), [], 1);  
        elseif ch == 2
            %G-channel
            result(G_Mos) = reshape(x(G_Mos), [], 1);
        elseif ch == 3
            %B-channel
            result(B_Mos) = reshape(x(B_Mos), [], 1); 
        end
        
    end    

return;
