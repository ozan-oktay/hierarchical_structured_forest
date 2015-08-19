function H = gradientHist3D(M, Gd, shrink, nOrients)

if (nOrients==0)
  H = [];
  return
end

if (nOrients~=6)
  error ('HoG implementation supports only 6 channels');
end

%load data.mat
%image = b / 5000;
%normRad = 4;       
%normConst = 0.05;
%grdSmooth = 3;
%shrink = 2;
% size of the input image
%[nrows, ncols, ndepth] = size(image);
% smoothing input image and gradient magnitude
%image_smoothed = convTri3D( image, grdSmooth );
%[M,~,Gd]       = gradientMag3D( image_smoothed, normRad, normConst ); clear image_smoothed;

% output, histogram bins and gradient direction normalization
[nrows, ncols, ndepth] = size(M);
bins= [1,0,0; -1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];
H   = zeros (nrows,ncols,ndepth,6,'single');
Gd  = Gd./repmat(sqrt(sum(Gd.^2,4)),[1,1,1,3]);

% Soft-binning - find the closest two bins
dotscores = zeros(nrows,ncols,ndepth,size(bins,1),'single');
for binId = 1 : size(bins,1)
   dotscores(:,:,:,binId) =     sum(   Gd(:,:,:,1) .* bins(binId,1) ...
                                     + Gd(:,:,:,2) .* bins(binId,2) ...
                                     + Gd(:,:,:,3) .* bins(binId,3), 4);
end
[binScore,binId] = sort(dotscores,4,'descend'); clear dotscores;

% Keep the two bins and discard the rest
binScore(:,:,:,3:end) = []; binScore = binScore ./ repmat(sum(binScore,4),[1,1,1,2]);
binId   (:,:,:,3:end) = [];

% Assign M values to the bins / soft orientation binning
for nn = 1:nrows
    for cc = 1:ncols
        for dd = 1:ndepth
            H(nn,cc,dd,binId(nn,cc,dd,1)) = binScore(nn,cc,dd,1) * M(nn,cc,dd);
            H(nn,cc,dd,binId(nn,cc,dd,2)) = binScore(nn,cc,dd,2) * M(nn,cc,dd);
        end
    end
end

% Spatial Pooling 
H = convTri3D( H, 2 ); 
H = imResample3D(H,size(H)/shrink);

end
