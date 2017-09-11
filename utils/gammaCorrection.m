function imgOut = gammaCorrection(imgIn, gammaTable)
%GAMMA Summary of this function goes here
%   Detailed explanation goes here
    [height, width] = size(imgIn);
    imgOut = zeros([height, width]);
    for y=1:height
        for x=1:width
            pos0 = floor(double(imgIn(y,x))/8) + 1 ;
            pos1 = floor(double(imgIn(y,x))/8) + 2 ;
            pos0 = max(min(pos0, 512), 1) ;
            pos1 = max(min(pos1, 512), 1) ;
            d0 = gammaTable(pos0);
            d1 = gammaTable(pos1);
            w1 = imgIn(y,x) + 8 - pos0*8 ;
            w0 = 8 - w1 ; 
            imgOut(y,x) = floor(double(d0*w0 + d1*w1)/2);
        end
    end 
    imgOut = uint16(imgOut);
end

