function [ I_conx ] = contextualFeatures3D( I_input, opts )
%CONTEXTUALFEATURES Summary of this function goes here
%   Detailed explanation goes here

if(opts.nChnsConx~=6)
  error('number of channels error (6)')
end

%Search neighborhood
r_patch = 1;
r_disp  = round(opts.imWidth/4);
xyz_disp = []; xyz_patch = [];

for radiusId = 1:numel(r_disp), r = r_disp(radiusId);
  xyz_t = [ r,0,0; -r,0,0; 0,r,0; 0,-r,0; 0,0,r; 0,0,-r ];
  xyz_disp = cat(1,xyz_disp,xyz_t);      
end

for r = 1:r_patch 
  xyz_t = [ r,0,0; 0,r,0; 0,0,r; -r,0,0; 0,-r,0; 0,0,-r];
  xyz_patch = cat(1,xyz_patch,xyz_t);      
end

%Feature binary string
size_input = size(I_input);
type_input = class(I_input);
I_conx     = zeros([size_input,size(xyz_disp,1)],type_input);
median_img = zeros([size_input,size(xyz_patch,1)]);

%Find the median image
for r = 1:size(xyz_patch,1)
  median_img(:,:,:,r) = volshift(I_input,xyz_patch(r,1),xyz_patch(r,2),xyz_patch(r,3));  
end
median_img = median(median_img,4);

%Pad 2D array
p          = [max(r_disp),max(r_disp),max(r_disp)];
median_img = padarray(median_img,p,'replicate','both');

%Loop over each xy_disp
for r = 1:size(xyz_disp,1)
  dep_m = (1:size_input(3))+max(r_disp);
  col_m = (1:size_input(2))+max(r_disp);
  row_m = (1:size_input(1))+max(r_disp);
  
  x_disp = xyz_disp(r,1);
  y_disp = xyz_disp(r,2);
  z_disp = xyz_disp(r,3);
  
  tmp_img         = median_img(row_m+y_disp,col_m+x_disp,dep_m+z_disp);
  I_conx(:,:,:,r) = cast(tmp_img - median_img(row_m,col_m,dep_m),type_input);
end

end

function im1shift=volshift(im1,x,y,z)

[m,n,o]=size(im1);

im1shift=im1;
x1s=max(1,x+1);
x2s=min(n,n+x);

y1s=max(1,y+1);
y2s=min(m,m+y);

z1s=max(1,z+1);
z2s=min(o,o+z);

x1=max(1,-x+1);
x2=min(n,n-x);

y1=max(1,-y+1);
y2=min(m,m-y);

z1=max(1,-z+1);
z2=min(o,o-z);

im1shift(y1:y2,x1:x2,z1:z2)=im1(y1s:y2s,x1s:x2s,z1s:z2s);
end
