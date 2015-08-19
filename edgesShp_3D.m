function  [chnsShp] = edgesShp_3D(pem,vtk,w2i,opts)
%EDGESSHP_3D Summary of this function goes here
%   Detailed explanation goes here

% IF INPUT IS AN EMPTY ARRAY THEN RETURN AN EMPTY ARRAY
if ( isempty(pem) || isempty(vtk) || isempty(w2i) ), chnsShp=[]; return; end
  
% MAP THEM TO IMAGE DOMAIN & FIND THE MEAN POINT 
vtk     = padarray(vtk,[0,1],1.0,'post');
vol_ind = w2i * vtk'; 
cnt_ind = round(mean(vol_ind,2));

% CROP THE WINDOW BASED ON THE SPECIFIED OPTS PARS
shpShrink = opts.shpShrink;
shpWidth  = opts.shpWidth;
shpDepth  = opts.shpDepth;
w         = opts.shpWidth/2;
d         = opts.shpDepth/2;
cnt       = cnt_ind + 1; % Matlab Indexing
pem       = convTri3D( pem, opts.shpSmooth );
pem       = pem( cnt(1)-w+1:cnt(1)+w, cnt(2)-w+1:cnt(2)+w, cnt(3)-d+1:cnt(3)+d );
chnsShp   = imResample3D(pem,size(pem)/shpShrink);
chnsShp   = single(chnsShp);
expSize   = [shpWidth/shpShrink,shpWidth/shpShrink,shpDepth/shpShrink];
assert (all(expSize==size(chnsShp)));

end

