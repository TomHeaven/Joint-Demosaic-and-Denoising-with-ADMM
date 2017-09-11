function [psnr, psnr_r, psnr_g, psnr_b] = comppsnr_rgb(I_hat, I)
% red channel
psnr_r = comppsnr(double(I_hat(:,:,1)), double(I(:,:,1)));
% green channel
psnr_g = comppsnr(double(I_hat(:,:,2)), double(I(:,:,2)));
% blue channel
psnr_b = comppsnr(double(I_hat(:,:,3)), double(I(:,:,3)));
psnr = (psnr_r + psnr_g + psnr_b) / 3;

end