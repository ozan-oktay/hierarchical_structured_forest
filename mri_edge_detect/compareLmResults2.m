function compareLmResults2()
%COMPARERESULTS Summary of this function goes here
%   Detailed explanation goes here
clear all; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath;
workdirectory  = '/vol/biomedic/users/oo2113/str_hier_forest_mri'; addpath(workdirectory);
patientMatFile = strcat(workdirectory,'/mritestingdata/dofs/patientParam.mat');
firstTxtFile   = strcat(workdirectory,'/mritestingdata/results/distances_new_prn016_shape_large.txt');
secondTxtFile  = strcat(workdirectory,'/mritestingdata/results/distances_prn000.txt');
addpath(genpath(horzcat(workdirectory,'/toolbox')));
  
firstFileID  = fopen(firstTxtFile, 'r');
secondFileID = fopen(secondTxtFile, 'r');
load(patientMatFile);

nLandmarks = 6;
formatSpec = ['%s %f %f %f %f %f %f %s %f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f'];
firstData  = textscan(firstFileID, formatSpec,'Delimiter',' ');
secondData = textscan(secondFileID, formatSpec,'Delimiter',' ');
errorsLm1  = cell(nLandmarks+1+nLandmarks*3,1);
errorsLm2  = cell(nLandmarks+1+nLandmarks*3,1);
rotations  = [];
scales     = [];

for firstPointInd=1:numel(firstData{1})
  firstPointLongName  = firstData{1}(firstPointInd);
  secondPointLongName = secondData{1}(firstPointInd);
  
  firstPointName=strsplit(firstPointLongName{1},'/');
  firstPointName=remSpace(firstPointName{end});
  secondPointName=strsplit(secondPointLongName{1},'/');
  secondPointName=remSpace(secondPointName{end});
  
  if (~strcmp(firstPointName,secondPointName)), error('the filenames are not matching\n'); end;
  
  count=0;
  for lmInd=1:numel(firstData)
    if(iscell(firstData{lmInd}(1))), continue; end; count = count+1;
    errorsLm1{count} = [errorsLm1{count}; firstData{lmInd}(firstPointInd)];
    errorsLm2{count} = [errorsLm2{count}; secondData{lmInd}(firstPointInd)];
  end
  
  for rotPointInd=1:numel(z_rot)
    foundflag=false;
    rotPointName=strsplit(filename(rotPointInd,:),'/');
    if (isequal(remSpace(rotPointName{end}),firstPointName))
      rotations = [rotations ; z_rot(rotPointInd)];
      scales    = [scales ; scale(rotPointInd)];
      foundflag = true;
    end
    if (foundflag), break; end
  end
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert the cell struct to an array
errorsLm1 = squeeze(cell2array(errorsLm1))-1;
errorsLm2 = squeeze(cell2array(errorsLm2))+1;

fprintf ('*/*/*/*/*/ Final Results /*/*/*/*/*/ \n');
fprintf (sprintf('\n*/*/*/* total number of matched images: %d */*/*/*\n',size(errorsLm1,1)));
fprintf ('distances txt name: %s\n',firstTxtFile);
for pointInd=1:nLandmarks
  fprintf (sprintf('*/*/*/ Landmark %d: mean: %6.3f std: %6.3f median: %6.3f \n',pointInd,mean(errorsLm1(:,pointInd)),std(errorsLm1(:,pointInd)),median(errorsLm1(:,pointInd))));
end
fprintf (sprintf('*/*/*/ Landmark C: mean: %6.3f std: %6.3f median: %6.3f \n',mean(errorsLm1(:,nLandmarks+1)),std(errorsLm1(:,nLandmarks+1)),median(errorsLm1(:,nLandmarks+1))));
fprintf ('================================================================\n');
fprintf (sprintf('*/*/*/ Landmark A: mean: %6.3f std: %6.3f median: %6.3f \n',mean(mean(errorsLm1(:,1:nLandmarks))),mean(std(errorsLm1(:,1:nLandmarks))),mean(median(errorsLm1(:,1:nLandmarks)))));
fprintf ('\n');
fclose(firstFileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show the distribution of the mean landmark errors wrt subjects
error_bins  = 0:1:30;
data        = mean(errorsLm1(:,1:nLandmarks),2);
pdfMean     = mean(data);
pdfStd      = std(data);
pdfRange    = 0:0.1:max(error_bins);
n_elem      = histc(data,error_bins);
norm_n_elem = n_elem / sum(n_elem);
c_elements  = cumsum(n_elem) / sum(n_elem);
figure(4); 
subplot(1,2,1);bar(error_bins,norm_n_elem,'BarWidth',1); xlim([0,max(error_bins)]); xlabel('Mean Landmark Error [mm]','FontSize',12); ylabel('Frequency','FontSize',12); grid on; hold on; plot(pdfRange,normpdf(pdfRange,pdfMean,pdfStd),'-r','Linewidth',3); 
h_legend=legend('Stratified Forest','Location','NorthEast'); set(h_legend,'FontSize',14);
title('Histogram of the Errors','FontSize',12);hold off;
subplot(1,2,2);plot(error_bins,c_elements,'-*b','LineWidth',2);  xlim([0,max(error_bins)]); xlabel('Mean Landmark Error [mm]','FontSize',12); ylabel('Percentage','FontSize',12); grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show the distribution of the mean landmark errors wrt subjects
error_bins  = 0:1:30;
data        = mean(errorsLm2(:,1:nLandmarks),2);
pdfMean     = mean(data);
pdfStd      = std(data);
pdfRange    = 0:0.1:max(error_bins);
n_elem      = histc(data,error_bins);
norm_n_elem = n_elem / sum(n_elem);
c_elements  = cumsum(n_elem) / sum(n_elem);
figure(4); hold on;
subplot(1,2,2);plot(error_bins,c_elements,'-*g','LineWidth',2);  xlim([0,max(error_bins)]); xlabel('Mean Landmark Error [mm]','FontSize',12); ylabel('Percentage','FontSize',12); grid on;
h_legend=legend('Stratified Forest','Hough Forest','Location','SouthEast'); set(h_legend,'FontSize',14);
title('Cumulative Distribution of the Errors','FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show the errors for each component ! 
fprintf ('================================================================\n');
fprintf (sprintf('*/*/*/ Landmark X: mean: %6.3f std: %6.3f median: %6.3f \n',mean(mean(errorsLm1(:,nLandmarks+1+1:3:end))),mean(std(errorsLm1(:,nLandmarks+1+1:3:end))),mean(median(errorsLm1(:,nLandmarks+1+1:3:end)))));
fprintf (sprintf('*/*/*/ Landmark Y: mean: %6.3f std: %6.3f median: %6.3f \n',mean(mean(errorsLm1(:,nLandmarks+1+2:3:end))),mean(std(errorsLm1(:,nLandmarks+1+2:3:end))),mean(median(errorsLm1(:,nLandmarks+1+2:3:end)))));
fprintf (sprintf('*/*/*/ Landmark Z: mean: %6.3f std: %6.3f median: %6.3f \n',mean(mean(errorsLm1(:,nLandmarks+1+3:3:end))),mean(std(errorsLm1(:,nLandmarks+1+3:3:end))),mean(median(errorsLm1(:,nLandmarks+1+3:3:end)))));
fprintf ('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function outstr = remSpace(str)
outstr= str(~isspace(str));
end

