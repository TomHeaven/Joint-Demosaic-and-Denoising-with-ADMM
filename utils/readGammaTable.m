function [ gammaTable ] = readGammaTable( gammaTableFileName )
%READGAMMATABLE Summary of this function goes here
%   gammaTableFileName = 
%   ['ipipe.txt_GammaR.txt' ; 'ipipe.txt_GammaG.txt'; 'ipipe.txt_GammaB.txt'];
%   512 => 1024

if ~exist('gammaTableFileName', 'var')
    gammaTableFileName = ['ipipe.txt_GammaR.txt' ; 'ipipe.txt_GammaG.txt'; 'ipipe.txt_GammaB.txt'];
end

gammaTable = zeros(3,512);
for t=1:3
    fid = fopen(gammaTableFileName(t,:), 'r');
    tmp = fscanf(fid, '%d', [2,512]);
    gammaTable(t, :) = tmp(1, :);
    fclose(fid);
end

end