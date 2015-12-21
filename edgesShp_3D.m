function  [chnsShp] = edgesShp_3D(pem,vtk,w2i,opts)
%EDGESSHP_3D Summary of this function goes here
%   Detailed explanation goes here

% IF INPUT IS AN EMPTY ARRAY THEN RETURN AN EMPTY ARRAY
if ( isempty(pem) || isempty(vtk) || isempty(w2i) ), chnsShp=[]; return; end
  
% MAP THEM TO IMAGE DOMAIN & FIND THE MEAN POINT 
vtk     = padarray(vtk,[0,1],1.0,'post');
vol_ind = w2i * vtk'; 
cnt_ind = round(mean(vol_ind,2));
chnsShp = cell(1,3);

%% COMPUTE THE LANDMARK DISTANCES AND DISTANCE RATIOS
num_landmarks    = size(vtk,1);
pairs            = nchoosek(1:num_landmarks,2);
num_combinations = size(pairs,1);
pair_distances   = zeros(num_combinations,1);

for loopId=1:num_combinations, cmb=pairs(loopId,:);
  pnt1 = vtk(cmb(1),:);
  pnt2 = vtk(cmb(2),:);
  pair_distances(loopId) = sqrt(sum( (pnt1(:) - pnt2(:)).^2 ));
end

num_distances    = size(pair_distances,1);
pairs            = nchoosek(1:num_distances,2);
num_combinations = size(pairs,1);
pair_ratios      = zeros(num_combinations,1);
for loopId=1:num_combinations, cmb=pairs(loopId,:);
  dist1 = pair_distances(cmb(1));
  dist2 = pair_distances(cmb(2));
  pair_ratios(loopId) = dist1 / dist2;
end
chnsShp{1} = single([pair_distances(:);pair_ratios(:)]);

%% CROP THE WINDOW BASED ON THE SPECIFIED OPTS PARS
shpShrink = opts.shpShrink;
shpWidth  = opts.shpWidth;
shpDepth  = opts.shpDepth;
w         = opts.shpWidth/2;
d         = opts.shpDepth/2;
cnt       = cnt_ind + 1; % Matlab Indexing
pem       = convTri3D( pem, opts.shpSmooth );
pem       = pem( cnt(1)-w+1:cnt(1)+w, cnt(2)-w+1:cnt(2)+w, cnt(3)-d+1:cnt(3)+d );

%% COMPUTE THE HoG FEATURES
fd_size  = (3*opts.nShpOrients+5) * floor(shpWidth/opts.nShpBinSize) * floor(shpWidth/opts.nShpBinSize);
fd_array = zeros(fd_size,shpDepth/shpShrink,'single');
for loopId = 1:(shpDepth/shpShrink)
  sliceId   = 1 + (loopId-1)*shpShrink;
  slice     = single(pem(:,:,sliceId)) / opts.ctmaxval;
  fd        = fhog( slice, opts.nShpBinSize, opts.nShpOrients, 0.2, false );
  %figure(1);im(slice);
  %figure(2);im(hogDraw(fd,25,1));
  fd_array(:,loopId) = fd(:);
end
chnsShp{2}  = single(fd_array(:));

%% Compute the appearance features
pemDown    = imResample3D(pem,size(pem)/shpShrink);
expSize    = [shpWidth/shpShrink,shpWidth/shpShrink,shpDepth/shpShrink]; assert (all(expSize==size(pemDown)));
chnsShp{3} = single(pemDown(:));
chnsShp    = [chnsShp{1};chnsShp{2};chnsShp{3}];

end
