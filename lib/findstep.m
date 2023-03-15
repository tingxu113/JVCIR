% This is a simple script to  automatically find the steps of the new priority definition
% by L. -J. Deng (UESTC)
% email: liangjian1987112@126.com
%% ==================================
function teststep = findstep(imgFilename, sourceRegion, w)
W = imread(imgFilename);
W = rgb2gray(double(W)/255);
BW2 = edge(W,'canny');  %prewitt; canny;
[m, n] = size(W);

a = nnz(BW2.*sourceRegion);   % non-zeros for out of mask
b =  m*n - nnz(~sourceRegion);   % total number for out of mask
% 
% a = nnz(BW2);
% b = m*n;

 
ration = w*a/b;  


nnz_maskregion = nnz(~sourceRegion);
%halfpatch = (2*w+1)*w;
halfpatch = 81/2;

est_totalstep = nnz_maskregion/halfpatch;

%teststep = ration*est_totalstep/(1+ration);

teststep = ration*est_totalstep;
teststep = floor(teststep);
% step.structure is our desired!!!
% from: ration = step.structure/(totalstep - step.structure)
