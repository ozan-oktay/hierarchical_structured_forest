function compareConfidenceRes()
%COMPARERESULTS Summary of this function goes here
%   Detailed explanation goes here

workdirectory = '/vol/biomedic/users/oo2113/str_hier_forest_mri'; addpath(workdirectory);
confTxtFile   = strcat(workdirectory,'/mritestingdata/results/confEstimations_new_prn016_shape_large.txt');
distTxtFile   = strcat(workdirectory,'/mritestingdata/results/distances_new_prn016_shape_large.txt');
addpath(genpath(horzcat(workdirectory,'/toolbox')));

confFileID     = fopen(confTxtFile, 'r');
confformatSpec = ['%s %s %s %f'];
confiData      = textscan(confFileID, confformatSpec,'Delimiter',' ');
valueInd       = [4];
confvalues     = cell(numel(valueInd),1);

nLandmarks     = 6;
distformatSpec = ['%s %f %f %f %f %f %f %s %f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f'];
distFileID     = fopen(distTxtFile, 'r');
distData       = textscan(distFileID, distformatSpec,'Delimiter',' ');
errorsLm       = cell(nLandmarks+1+nLandmarks*3,1);

for firstPointInd=1:numel(confiData{2})
    
  for lmInd=1:size(confvalues,1), iter=valueInd(lmInd);
    confvalues{lmInd} = [confvalues{lmInd}; confiData{iter}(firstPointInd)];
  end
  
  fistPointLongName=confiData{2}(firstPointInd);
  firstPointName=strsplit(fistPointLongName{1},'/');
  firstPointName=remSpace(firstPointName{end});
  
  for secondPointInd=1:numel(distData{1})
    foundflag=false;
    distPointFullName=distData{1}(secondPointInd,:);
    distPointName=strsplit(distPointFullName{1},'/');
    if (isequal(remSpace(distPointName{end}),firstPointName)),foundflag = true; end
    if (foundflag), break; end
  end
  
  count=0;
  for lmInd=1:numel(distData)
    if(iscell(distData{lmInd}(1))), continue; end; count = count+1;
    errorsLm{count} = [errorsLm{count}; distData{lmInd}(secondPointInd)];
  end
end

errorsLm      = squeeze(cell2array(errorsLm)) - 1;
confvalues    = cell2array(confvalues);
meanErrorPSub = mean(errorsLm(:,1:nLandmarks),2);
[RHO,PVAL] = corr(confvalues,meanErrorPSub,'type','Pearson');

figure(1);
plot(confvalues,meanErrorPSub,'*b');
xlabel('Confidence Values');
ylabel('Mean Landmark Error');
title(sprintf('Rho Val %f, P Val %f',RHO,PVAL));
grid on;

end


function outstr = remSpace(str)
outstr= str(~isspace(str));
end