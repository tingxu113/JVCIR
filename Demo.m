% This is a demo for "Exemplar-based image inpainting using adaptive 
% two-stage structure-tensor based priority function and nonlocal filtering" 
% Journal of Visual Communication and Image Representation, 2022.
% 
% Created by Ting Xu (UESTC)
% Updated: Mar. 15, 2023.
% ==================================================
% If you find this code helpful, please kindly cite this article.
% ==================================================
%% ============= set path===============
clear;close all
addpath('TestImages');
addpath('lib');
ww = 'angle';
dirTe = strcat('TestImages\', ww, '1.png');
dirTe2 = strcat('TestImages\', ww, '.png');
color = [255 0 0];
%% Main function
[i1,i2,i3, c,d,t,mov,conf,deta]=inpaintNewDefi_ting(dirTe2,dirTe,color);
%% ============= present results============%
figure;subplot(131);imshow(uint8(i2)); title('Original image');
       subplot(132);imshow(uint8(i3)); title('Fill region');
       subplot(133);imshow(uint8(i1)); title('Inpainted image');

 for i=1:3
     SSIMVAL(i) = ssim(i1(:,:,i), i2(:,:,i));
 end
     SSIMVAL=1/3*sum(SSIMVAL(1:3))
for i = 1:3
    X1 = i1(:,:,i);
    i1(:,:,i)=X1./max(max(X1(:)));
end
for i = 1:3
    X2 = i2(:,:,i);
    i2(:,:,i)=X2./max(max(X2(:)));
end
psnr_p=psnr3(i1/255,i2/255)
disp('------------------------------------------------')
display(['Total time','    ',num2str(t(1)+t(2)+t(3)+t(4)+t(5))])
disp('--------------   RESULTS TABLE   --------------')
disp('Ready phase  | Compute gradient  |Get the priorities  |  Get the exemplars  | Update and copy ')
display([num2str(t(4)), '                ',num2str(t(5)), '               ',num2str(t(1)), '                      ', num2str(t(2)), '                      ', num2str(t(3))])




