% This is a main function for "Exemplar-based image inpainting using adaptive 
% two-stage structure-tensor based priority function and nonlocal filtering" 
% Journal of Visual Communication and Image Representation, 2022.
% 
% Created by Ting Xu (UESTC)
% Updated: Mar. 15, 2023.
% ==================================================
% If you find this code helpful, please kindly cite this article.
%% ===========================================
function [inpaintedImg,origImg,fillImg,C,D,t,ConfMovie,DetaMovie,ratio] = inpaintNewDefi_ting(imgFilename,fillFilename,fillColor)
%INPAINT  Exemplar-based inpainting.
% Usage:   [inpaintedImg,origImg,fillImg,C,D,fillMovie] ...
%                = inpaint(imgFilename,fillFilename,fillColor)
% Inputs: 
%   imgFilename    Filename of the original image.
%   fillFilename   Filename of the image specifying the fill region. 
%   fillColor      1x3 RGB vector specifying the color used to specify
%                  the fill region.
%   opts           parameter
% Outputs:
%   inpaintedImg   The inpainted image; an MxNx3 matrix of doubles. 
%   origImg        The original image; an MxNx3 matrix of doubles.
%   fillImg        The fill region image; an MxNx3 matrix of doubles.
%   C              MxN matrix of confidence values accumulated over all iterations.
%   D              MxN matrix of data term values accumulated over all iterations.
%   fillMovie      A Matlab movie struct depicting the fill region over
%   time. 
%------------------------------------------------------------------------
t1=0; t2=0; t3=0; t4=0; t5=0; 
tic
[img,fillImg,fillRegion] = loadimgs(imgFilename,fillFilename,fillColor);
img = double(img);
origImg = img;
ind = img2ind(img);
sz = [size(img,1) size(img,2)];
sourceRegion = ~fillRegion;

%mask's ratio in image
ratio= sum(sum(fillRegion))/(sz(1)*sz(2));
s=ratio*sz(1)*sz(2);
teststep=findstep(imgFilename, sourceRegion,17);%17 can be adjusted

% Initialize isophote values
[Ix(:,:,3), Iy(:,:,3)] = gradient(img(:,:,3));
[Ix(:,:,2), Iy(:,:,2)] = gradient(img(:,:,2));
[Ix(:,:,1), Iy(:,:,1)] = gradient(img(:,:,1));
Ix = sum(Ix,3)/(3); Iy = sum(Iy,3)/(3);
   
% Initialize confidence and data terms
C = double(sourceRegion);
D = repmat(-.1,sz);
iter = 1;

t4=t4+toc;
% Visualization stuff
if nargout==9
fillMovie(1).cdata=uint8(img); 
fillMovie(1).colormap=[]; 
E(:,:,1)=C;E(:,:,2)=C;E(:,:,3)=C;
ConfMovie(1).cdata=uint8(255*E);
ConfMovie(1).colormap=[];
F(:,:,1)=D;F(:,:,2)=D;F(:,:,3)=D;
DetaMovie(1).cdata=uint8(255*F);
DetaMovie(1).colormap=[];
origImg(1,1,:) = fillColor;
end

rand('state',0);

    % Loop until entire fill region has been covered
while any(fillRegion(:))  
     tic
    % Find contour & normalized gradients of fill region
    fillRegionD = double(fillRegion); % Marcel 11/30/05
    dR = find(conv2(fillRegionD,[1,1,1;1,-8,1;1,1,1],'same')>0); % Marcel 11/30/05
    [Nx,Ny] = gradient(double(~fillRegion)); % Marcel 11/30/05
    N = [Nx(dR(:)) Ny(dR(:))];    
    N = normr(N);  
    N(~isfinite(N))=0; % handle NaN and Inf

    t5 = t5 +toc;
    tic
for k=dR';  Hp = getpatch(sz,k); q = Hp(~(fillRegion(Hp))); C(k) = sum(C(q))/numel(Hp);  
      %Compute data term 
    c1=ones(size(q,1),1);
    r1=ones(size(q,1),1);
    for i=1:size(q,1)
              x1=q(i)/sz(1);
              if x1==round(x1)
              c1(i)=x1;r1(i)=sz(1);
              else
              c1(i)=floor(x1)+1;r1(i)=q(i)-sz(1)*floor(x1);
              end
    end
    bb=size(q,1);
    ww=zeros(bb,1);
    sigma=0.3;%adjust for better performance
        
    for j=1:size(q,1)
        for kk=1:size(q,1)
              s(j)=(1/2*pi*sigma^2)*(exp(-((c1(kk)-c1(j))^2+(r1(kk)-r1(j))^2)/2*sigma^2));
              ww(j)=ww(j)+s(j);
        end
        ww(j)=ww(j)-1/2*pi*sigma^2;
    end   
    G=sum((ww.*Ix(q).*Ix(q)+ww.*Iy(q).*Iy(q))/max(ww));
    D(k)=G*norm(N,2);
end
      
if (iter<teststep)
  priorities = D(dR);
else
  priorities = C(dR);
end

    t1 = t1 +toc;

    tic
      [~,ndx] = max(priorities(:));
      p = dR(ndx(1));
      [Hp,rows,cols] = getpatch(sz,p);
      toFill = fillRegion(Hp);
      inputImg = img;   
 kkk=2;
for i = 1:kkk
[Ep,erows,ecols,w2] = getbigpatch(sz,p);
sRegion = sourceRegion(Ep)';  
bigpatchimg = inputImg(erows,ecols,:); 
[Hq,Hq1,best] = newbestexemplar(inputImg,bigpatchimg,img(rows,cols,:),toFill',sRegion,sz,p,w2);
cadpat{i}=Hq1;
cadpat1{i}=Hq;
%location{i}=best;
inputImg(best(1)+15:best(2),(best(3)+15:best(4))',:) = 0;
end

     t2 = t2 +toc;
     tic
     toFill = logical(toFill);  % Marcel 11/30/05
     fillRegion(Hp(toFill)) = false;

      C(Hp(toFill))  = C(p);

      fil = repmat(~toFill', [1 1 3]);
      gama2=0.1;%adjust for better performance
      s1=[];
      s2=fil.*(origImg(rows,cols,:)./255);

      for i=1:kkk
        s1{i}=fil .*(cadpat{i}./255);
       b(i)=(norm(s1{i}(:)-s2(:),2))^2;
       o(i)=exp(-b(i)/2*gama2);
      end  
       s=sum(o(1:kkk)); 

      Hp1=[];
      for i=1:kkk
          ind(Hp(toFill)) = ind(cadpat1{i}(toFill));  
          img(rows,cols,:) = ind2img(ind(rows,cols),origImg);
          Hp1{i}=img(rows,cols,:);
      end
      img(rows,cols,:) =o(1)/s* Hp1{1}+o(2)/s* Hp1{2};


    t3 = t3 +toc;
      if nargout==9
        ind2 = ind;
        ind2(logical(fillRegion)) = 1;          % Marcel 11/30/05

        fillMovie(iter).cdata=uint8(ind2img(ind2,origImg)); 
        fillMovie(iter).colormap=[];
        E(:,:,1)=C;E(:,:,2)=C;E(:,:,3)=C;
        ConfMovie(iter).cdata=uint8(255*E);
        ConfMovie(iter).colormap=[];
        F(:,:,1)=D;F(:,:,2)=D;F(:,:,3)=D;
        DetaMovie(iter).cdata=uint8(255*F);
        DetaMovie(iter).colormap=[];
      end

      iter = iter+1;
end
inpaintedImg=img;

t=[t1;t2;t3;t4;t5];


   
   







 
        


%---------------------------------------------------------------------
% Scans over the entire image (with a sliding window)
% for the exemplar with the lowest error. Calls a MEX function.
%---------------------------------------------------------------------
%  to get the new exemplar with patch-in-patch method

function [Hq,Hq1,best] = newbestexemplar(inputImg,img,Ip,toFill,sourceRegion,sz,p,w)
m=size(Ip,1); mm=size(img,1); n=size(Ip,2); nn=size(img,2);
best = bestexemplarhelper(mm,nn,m,n,img,Ip,toFill,sourceRegion);
g(2) = fix(p/sz(1))+1; g(1) = mod(p,sz(1));
ff=16;% note that ff must be the 2*w.
 if (g(1)<w+1)&&(g(2)<w+1)
     best(1) = best(1); best(2) = best(1)+ff;
     best(3) = best(3); best(4) = best(3)+ff;
 elseif (g(1)<w+1)&&(w+1<=g(2))
     best(1) = best(1); best(2) = best(1)+ff;
     best(3) = g(2)+best(3)-(w+1);  best(4) = best(3)+ff;
 elseif (w+1<=g(1))&&(g(2)<w+1)
     best(1) = g(1)+best(1)-(w+1);  best(2) = best(1)+ff;
     best(3) = best(3); best(4) = best(3)+ff;
 else
     best(1) = g(1)+best(1)-(w+1);  best(2) = best(1)+ff;
     best(3) = g(2)+best(3)-(w+1);  best(4) = best(3)+ff;
 end
Hq1 =inputImg(best(1):best(2),best(3):best(4),:); 
Hq = sub2ndx(best(1):best(2),(best(3):best(4))',sz(1));

function [Hp,rows,cols] =  getpatch(sz,p)
% [x,y] = ind2sub(sz,p);  % 2*w+1 == the patch size
w=8;%adjust for better performance
p=p-1; y=floor(p/sz(1))+1; p=rem(p,sz(1)); x=floor(p)+1;
rows = max(x-w,1):min(x+w,sz(1));
cols = (max(y-w,1):min(y+w,sz(2)))';
Hp = sub2ndx(rows,cols,sz(1));

% patch-in-patch�� to get bigger patch
function [Ep,erows,ecols,w] = getbigpatch(sz,p)
% [x,y] = ind2sub(sz,p);  % 2*w+1 == the patch size
w =81;% big patchsize
p=p-1; y=floor(p/sz(1))+1; p=rem(p,sz(1)); x=floor(p)+1;
erows = max(x-w,1):min(x+w,sz(1));
ecols = (max(y-w,1):min(y+w,sz(2)))';
Ep = sub2ndx(erows,ecols,sz(1));
%-----------------3------------------------------------------------
% Converts the (rows,cols) subscript-style indices to Matlab index-style
% indices.  Unforunately, 'sub2ind' cannot be used for this.
%---------------------------------------------------------------------
function N = sub2ndx(rows,cols,nTotalRows)
X = rows(ones(length(cols),1),:);
Y = cols(:,ones(1,length(rows)));
N = X+(Y-1)*nTotalRows;


%---------------------------------------------------------------------
% Converts an indexed image into an RGB image, using 'img' as a colormap
%---------------------------------------------------------------------
function img2 = ind2img(ind,img)
for i=3:-1:1, temp=img(:,:,i); img2(:,:,i)=temp(ind); end


%---------------------------------------------------------------------
% Converts an RGB image into a indexed image, using the image itself as
% the colormap.
%---------------------------------------------------------------------
function ind = img2ind(img)
s=size(img); ind=reshape(1:s(1)*s(2),s(1),s(2));


%---------------------------------------------------------------------
% Loads the an image and it's fill region, using 'fillColor' as a marker
% value for knowing which pixels are to be filled.
%---------------------------------------------------------------------
function [img,fillImg,fillRegion] = loadimgs(imgFilename,fillFilename,fillColor)
img = imread(imgFilename); fillImg = imread(fillFilename);
fillRegion = fillImg(:,:,1)==fillColor(1) & ...
    fillImg(:,:,2)==fillColor(2) & fillImg(:,:,3)==fillColor(3);

