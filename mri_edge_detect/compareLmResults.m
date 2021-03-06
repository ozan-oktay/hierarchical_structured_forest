function compareLmResults()
%COMPARERESULTS Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath;
workdirectory  = '/vol/medic02/users/oo2113/str_hier_forest_mri'; addpath(workdirectory);
patientMatFile = strcat(workdirectory,'/mristacomdata/dofs/patientParam.mat');
firstTxtFile   = strcat(workdirectory,'/mristacomdata/results/distances_interUser.txt');
addpath(genpath(horzcat(workdirectory,'/toolbox')));
  
firstFileID  = fopen(firstTxtFile, 'r');
load(patientMatFile);

nLandmarks = 6;
formatSpec = ['%s %f %f %f %f %f %f %s %f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f'];
firstData  = textscan(firstFileID, formatSpec,'Delimiter',' ');
errorsLm   = cell(nLandmarks+1+nLandmarks*3,1);
rotations  = [];
scales     = [];

for firstPointInd=1:numel(firstData{1})
  fistPointLongName=firstData{1}(firstPointInd);
  firstPointName=strsplit(fistPointLongName{1},'/');
  firstPointName=remSpace(firstPointName{end});
  
  count=0;
  for lmInd=1:numel(firstData)
    if(iscell(firstData{lmInd}(1))), continue; end; count = count+1;
    errorsLm{count} = [errorsLm{count}; firstData{lmInd}(firstPointInd)];
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
errorsLm           = squeeze(cell2array(errorsLm)) - 1;
[XY_set,XZ_set]    = computeErrorsInPlane(errorsLm,nLandmarks);
[InsertPnt_XY_set] = computeInsertPntErrorsInPlane(errorsLm,nLandmarks);

fprintf ('*/*/*/*/*/ Final Results /*/*/*/*/*/ \n');
fprintf (sprintf('\n*/*/*/* total number of matched images: %d */*/*/*\n',size(errorsLm,1)));
fprintf ('distances txt name: %s\n',firstTxtFile);
for pointInd=1:nLandmarks
  fprintf (sprintf('*/*/*/ Landmark %d: mean: %6.3f std: %6.3f median: %6.3f \n',pointInd,mean(errorsLm(:,pointInd)),std(errorsLm(:,pointInd)),median(errorsLm(:,pointInd))));
end
fprintf (sprintf('*/*/*/ Landmark C: mean: %6.3f std: %6.3f median: %6.3f \n',mean(errorsLm(:,nLandmarks+1)),std(errorsLm(:,nLandmarks+1)),median(errorsLm(:,nLandmarks+1))));
fprintf ('================================================================\n');
fprintf (sprintf('*/*/*/ Landmark A: mean: %6.3f std: %6.3f median: %6.3f \n',mean(mean(errorsLm(:,1:nLandmarks))),mean(std(errorsLm(:,1:nLandmarks))),mean(median(errorsLm(:,1:nLandmarks)))));
fprintf (sprintf('*/*/*/ Landmark A: mean: %6.3f std: %6.3f median: %6.3f (XY)\n',XY_set(1),XY_set(2),XY_set(3)));
fprintf (sprintf('*/*/*/ Landmark A: mean: %6.3f std: %6.3f median: %6.3f (XZ)\n',XZ_set(1),XZ_set(2),XZ_set(3)));
fprintf (sprintf('*/*/*/ Landmark A: mean: %6.3f std: %6.3f median: %6.3f (XY_Insert_Pnts)\n',InsertPnt_XY_set(1),InsertPnt_XY_set(2),InsertPnt_XY_set(3)));
fprintf ('\n');
fclose(firstFileID);

% plot the subplots of rotations vs landmark error
figure(1); 
for ind=1:nLandmarks, subplot(2,3,ind);   
  plot(rotations,errorsLm(:,ind),'*b'); [RHO,PVAL] = corr(rotations,errorsLm(:,ind),'type','Pearson'); rangeY=[0:0.1:35]; rangeX=[min(rotations),max(rotations)];
  meanY=mean(errorsLm(:,ind)); stdY=std(errorsLm(:,ind)); maxX=max(rotations); pdfY=normpdf(rangeY,meanY,stdY); hold on; plot(maxX-300*pdfY,rangeY,'-r','Linewidth',3); plot(maxX-300*max(pdfY),meanY,'*r','Linewidth',7); hold off;
  title(sprintf('rho %f, pval %f',RHO,PVAL)); grid on; ylabel('Landmark Errors in mm'); xlabel('Rotation around the z-axis in degrees'); ylim([rangeY(1) rangeY(end)]); xlim([rangeX(1) rangeX(2)])
end;
% plot the subplots of scales vs landmark error
figure(2); 
for ind=1:nLandmarks, subplot(2,3,ind);   
  plot(scales,errorsLm(:,ind),'*b'); [RHO,PVAL] = corr(scales,errorsLm(:,ind),'type','Pearson'); rangeY=[0:0.1:35]; rangeX=[min(scales),max(scales)];
  meanY=mean(errorsLm(:,ind)); stdY=std(errorsLm(:,ind)); maxX=max(scales); pdfY=normpdf(rangeY,meanY,stdY); hold on; plot(maxX-3*pdfY,rangeY,'-r','Linewidth',3); plot(maxX-3*max(pdfY),meanY,'*r','Linewidth',7); hold off;
  title(sprintf('rho %f, pval %f',RHO,PVAL)); grid on; ylabel('Landmark Errors in mm'); xlabel('Rotation around the z-axis in degrees'); ylim([rangeY(1) rangeY(end)]); xlim([rangeX(1) rangeX(2)])
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show the distribution of the mean landmark errors wrt subjects
error_bins  = 0:1:30;
data        = mean(errorsLm(:,1:nLandmarks),2);
pdfMean     = mean(data);
pdfStd      = std(data);
pdfRange    = 0:0.1:max(error_bins);
n_elem      = histc(data,error_bins);
norm_n_elem = n_elem / sum(n_elem);
c_elements  = cumsum(n_elem) / sum(n_elem);
figure(3); 
subplot(1,2,1);bar(error_bins,norm_n_elem,'BarWidth',1); xlim([0,max(error_bins)]); xlabel('Mean Landmark Error [mm]','FontSize',12); ylabel('Frequency','FontSize',12); grid on; hold on; plot(pdfRange,normpdf(pdfRange,pdfMean,pdfStd),'-r','Linewidth',3); hold off;
subplot(1,2,2);plot(error_bins,c_elements,'-*b','LineWidth',2);  xlim([0,max(error_bins)]); xlabel('Mean Landmark Error [mm]','FontSize',12); ylabel('Percentage','FontSize',12); grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show the errors for each component ! 
fprintf ('================================================================\n');
for pointInd=((1:nLandmarks*3)+nLandmarks+1)
  fprintf (sprintf('*/*/*/ Landmark %d: mean: %6.3f std: %6.3f median: %6.3f \n',pointInd,mean(errorsLm(:,pointInd)),std(errorsLm(:,pointInd)),median(errorsLm(:,pointInd))));
end
fprintf ('================================================================\n');
fprintf (sprintf('*/*/*/ Landmark X: mean: %6.3f std: %6.3f median: %6.3f \n',mean(mean(errorsLm(:,nLandmarks+1+1:3:end))),mean(std(errorsLm(:,nLandmarks+1+1:3:end))),mean(median(errorsLm(:,nLandmarks+1+1:3:end)))));
fprintf (sprintf('*/*/*/ Landmark Y: mean: %6.3f std: %6.3f median: %6.3f \n',mean(mean(errorsLm(:,nLandmarks+1+2:3:end))),mean(std(errorsLm(:,nLandmarks+1+2:3:end))),mean(median(errorsLm(:,nLandmarks+1+2:3:end)))));
fprintf (sprintf('*/*/*/ Landmark Z: mean: %6.3f std: %6.3f median: %6.3f \n',mean(mean(errorsLm(:,nLandmarks+1+3:3:end))),mean(std(errorsLm(:,nLandmarks+1+3:3:end))),mean(median(errorsLm(:,nLandmarks+1+3:3:end)))));
fprintf ('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function outstr = remSpace(str)
outstr= str(~isspace(str));
end

function [XY_set,XZ_set]=computeErrorsInPlane(errorsLm,nLandmarks)

XY_set  = zeros(3);
XZ_set  = zeros(3);

Xerrors = errorsLm(:,nLandmarks+1 + (1:3:nLandmarks*3) );
Yerrors = errorsLm(:,nLandmarks+1 + (2:3:nLandmarks*3) );
Zerrors = errorsLm(:,nLandmarks+1 + (3:3:nLandmarks*3) );

XYerrors = sqrt(Xerrors.^2 + Yerrors.^2);
XZerrors = sqrt(Xerrors.^2 + Zerrors.^2);

XY_set(1) = mean(mean(XYerrors));
XY_set(2) = mean(std(XYerrors));
XY_set(3) = mean(median(XYerrors));

XZ_set(1) = mean(mean(XZerrors));
XZ_set(2) = mean(std(XZerrors));
XZ_set(3) = mean(median(XZerrors));

end

function [XY_set]=computeInsertPntErrorsInPlane(errorsLm,nLandmarks)

XY_set  = zeros(3);
Xerrors = errorsLm(:,nLandmarks+1 + [1,7]);
Yerrors = errorsLm(:,nLandmarks+1 + [2,8]);

XYerrors = sqrt(Xerrors.^2 + Yerrors.^2);

XY_set(1) = mean(mean(XYerrors));
XY_set(2) = mean(std(XYerrors));
XY_set(3) = mean(median(XYerrors));

end
