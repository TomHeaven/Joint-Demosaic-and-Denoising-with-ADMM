% Copyright 2014 NVIDIA Corporation.  All Rights Reserved
% This software is for provided for non-commercial and research purposes only
% Paper: FlexISP: A Flexible Camera Image Processing Framework, SIGGRAPH
% Asia 2014
% WebSite: https://research.nvidia.com/publication/flexisp-flexible-camera-image-processing-framework

function [ result ] = observationMat( x, channels, flag  )
    %Do the image formation
    sizeX = size(x);
    sizeI = sizeX(1:2);
    
    result = zeros( sizeX );
    for ch = 1:size(x,3)
        
        %Iterate over the channels
        if flag > 0      

            % A
            result(:,:,ch) = demosaickingMatrix(  imfilter(x(:,:,ch), channels(ch).K, 'same', 'replicate'), ch, sizeI, flag  );                               

        elseif flag < 0 

            % A'
            result(:,:,ch) = imfilter( reshape( demosaickingMatrix( x(:,:,ch), ch, sizeI, flag  ), sizeI) , fliplr(flipud(channels(ch).K)), 'same', 'conv', 'replicate');  

        else

            % A'A
            Ax = demosaickingMatrix( imfilter(x(:,:,ch), channels(ch).K, 'same', 'replicate'), ch, sizeI, flag );                                         
            result(:,:,ch) = imfilter(Ax, fliplr(flipud(channels(ch).K)), 'same', 'conv', 'replicate');  

        end  
        
    end
    
return;