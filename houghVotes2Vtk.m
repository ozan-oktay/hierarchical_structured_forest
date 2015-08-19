% Written by Ozan Oktay 2015 
% Inputs
%   HV: Accumulated Hough Votes of size H x W x D x CH
% varargin
%   method:       maximum, meanshift 
%   savename:     output vtk filename in world coordinates 
%   refimagename: input reference imagename used for image2world mapping
%   

function houghVotes2Vtk (HV,varargin)

dfs={'method','max', 'savename',[], 'refimagename',[], 'spacing',[1,1,1]};
[method, savename, refimagename, spacing] = getPrmDflt(varargin,dfs,1);

method=find(strcmpi(method,{'max','meanshift','weighted'}))-1;
if(isempty(method)), error('unknown point selection criteria: %s',method); end
if(isempty(savename)||isempty(refimagename)), return; end;
if(numel(size(HV))<4), error('houghVotes_2_Vtk supports only 4D Hough images'); end;

% Collect information about the hough vote image
siz        = size(HV);
nLandmarks = siz(end);
nOutputs   = nLandmarks * (numel(siz)-1);
outputArr  = [];

%Select the maximum probability point in each channel  
if (method == 0)
  for ii=1:nLandmarks 
    hChannel  = HV(:,:,:,ii);
    [~,index] = max(hChannel(:));
    [r,c,s]   = ind2sub(siz(1:end-1),index);
    outputArr = single([outputArr; [r-1, c-1, s-1]]); 
  end
  
%Apply meanshift density fitting and select the maximum point based on that
elseif (method == 1)
    for ii=1:nLandmarks  
      % Normalize the Hough Votes and collect samples 
      hChannel = HV(:,:,:,ii); hChannel=hChannel/max(hChannel(:));
      hMask    = hChannel > 0.25;
      weights  = double (hChannel(hMask));
      indices  = find (hMask==true);
      [X,Y,Z]  = ind2sub (siz(1:3),indices);
      ms_data  = cat(2,X,Y,Z);
      ms_data  = ms_data.*repmat(spacing(1:3),size(ms_data,1),1); %size in mm
      scale    = sqrt(max(var(ms_data)));
      ms_data  = ms_data / scale;   

      % Mean-shift parameters
      radius   = 10 / scale;     % size in mm
      maxIter  = 50;             % maximum number of iterations
      rate     = 0.20;           % gradient ascent learning rate
      nThreads = 64;             % Number of threads used in computation

      % Run the mean-shift algorithm on discretized data 
      [~,M,~]   = meanShiftWeighted(ms_data,radius,weights,rate,maxIter,0,nThreads);
      means     = single((M .* scale) ./ repmat(spacing(1:3),size(M,1),1));
      outputArr = [outputArr; [means(1,1)-1,means(1,2)-1,means(1,3)-1]]; %#ok<AGROW>
    end
    
%Apply weighted averaging 
elseif (method == 2)
  for ii=1:nLandmarks 
    hChannel      = HV(:,:,:,ii);
    [r_g,c_g,s_g] = ndgrid(1:size(hChannel,1),1:size(hChannel,2),1:size(hChannel,3));
    r             = sum(hChannel(:).*r_g(:))/sum(hChannel(:));
    c             = sum(hChannel(:).*c_g(:))/sum(hChannel(:));
    s             = sum(hChannel(:).*s_g(:))/sum(hChannel(:));
    outputArr     = single([outputArr; [r-1, c-1, s-1]]); 
  end  
end

% Write the results to a txt and vtk file 
assert(numel(outputArr)==nOutputs);
pnt2Vtk(refimagename,savename,outputArr);

end
