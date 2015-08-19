function [ I_conx ] = contextualFeatures( I_input )
%CONTEXTUALFEATURES Summary of this function goes here
%   Detailed explanation goes here

%Search neighborhood
r_patch = 1;
r_disp  = [8,16];
xy_disp = []; xy_patch = [];

for radiusId = 1:numel(r_disp), r = r_disp(radiusId);
  xy_t = [ r,r; -r,r; -r,-r; r,-r ];
  xy_disp = cat(1,xy_disp,xy_t);      
end

for r = 1:r_patch 
  xy_t = [ r,0; r,r; 0,r; -r,r; ...
          -r,0; -r,-r; 0,-r; r,-r];
  xy_patch = cat(1,xy_patch,xy_t);      
end

%Feature binary string
size_input = size(I_input);
type_input = class(I_input);
I_conx     = zeros([size_input,size(xy_disp,1)],type_input);
median_img = zeros([size(xy_patch,1),size_input]);

%Find the median image
for r = 1:size(xy_patch,1)
  median_img(r,:,:) = imshift(I_input,xy_patch(r,1),xy_patch(r,2));  
end
median_img = permute(median(median_img,1),[2,3,1]);

%Pad 2D array
p             = [max(r_disp),max(r_disp)];
median_padded = padarray(median_img,p,'replicate','both');

%Loop over each xy_disp
for r = 1:size(xy_disp,1)
  
  col_m = (1:size_input(2))+max(r_disp);
  row_m = (1:size_input(1))+max(r_disp);
  
  x_disp = xy_disp(r,1);
  y_disp = xy_disp(r,2);
  
  tmp_img       = median_padded(row_m+y_disp,col_m+x_disp);
  I_conx(:,:,r) = cast(tmp_img - I_input,type_input);
end

end

function im1shift=imshift(im1,x,y)

[m,n,o]=size(im1);

im1shift=im1;
x1s=max(1,x+1);
x2s=min(n,n+x);

y1s=max(1,y+1);
y2s=min(m,m+y);

x1=max(1,-x+1);
x2=min(n,n-x);

y1=max(1,-y+1);
y2=min(m,m-y);

im1shift(y1:y2,x1:x2,:)=im1(y1s:y2s,x1s:x2s,:);
end