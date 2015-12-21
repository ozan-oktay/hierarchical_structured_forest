function comparePoseResults()
%COMPARERESULTS Summary of this function goes here
%   Detailed explanation goes here

firstTxtFile    = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mristacomdata/results/poseEstimations_new_prn016_shape_large.txt';
workdirectory   = '/vol/biomedic/users/oo2113/str_hier_forest_mri'; addpath(workdirectory);
addpath(genpath(horzcat(workdirectory,'/toolbox')));
addpath(genpath('/homes/oo2113/workspace/matlab/plotTools/violin'));

firstFileID  = fopen(firstTxtFile, 'r');

formatSpec= ['%s %s %s %f %f %s %f %f %s %f %f'];
firstData = textscan(firstFileID, formatSpec,'Delimiter',' ');
valueInd  = [4,5,7,8,10,11];
values    = cell(numel(valueInd),1);

for firstPointInd=1:numel(firstData{1})
  for lmInd=1:size(values,1), iter=valueInd(lmInd);
    values{lmInd} = [values{lmInd}; firstData{iter}(firstPointInd)];
  end
end

% mean error for pose 
meanRotErr    = mean( sqrt( (values{1} - values{3}).^2 ) );
meanScaleErr  = mean( sqrt( (values{2} - values{4}).^2 ) );
meanRotConf   = mean(values{5});
meanScaleConf = mean(values{6});
stdRotErr     = std( sqrt( (values{1} - values{3}).^2) );
stdScaleErr   = std( sqrt( (values{2} - values{4}).^2) );

fprintf ('*/*/*/*/*/ Final Results  \n');
fprintf (sprintf('*/*/*/*/*/ Total number of images: %d \n',size(values{1},1)));
fprintf ('*/*/*/*/*/ Txt name: %s\n',firstTxtFile);
fprintf ('================================================================================================================================\n');
rotationText = sprintf('/*/*/*/*/* Mean Rotation Err: %3.4f -- Std: %3.4f -- Mean Conf Val: %3.4f */*/*/*/*/\n',meanRotErr,stdRotErr,meanRotConf);     fprintf(rotationText);
scaleText    = sprintf('/*/*/*/*/* Mean Scale    Err: %3.4f -- Std: %3.4f -- Mean Conf Val: %3.4f */*/*/*/*/\n',meanScaleErr,stdScaleErr,meanScaleConf); fprintf(scaleText);
fprintf ('\n');
fclose(firstFileID);

% display the distribution of the rotation and scale values in testing set
close all;
grndRot   = values{1}; figure(1); violin(grndRot); ylabel('Rotation Values','FontSize',14); grid on; title(rotationText,'FontSize',10); 
grndScale = values{2}; figure(2); violin(grndScale); ylabel('Scale Values','FontSize',14); grid on;  title(scaleText,'FontSize',10);    

% Range of grndRotation 
minGrndRot = min(grndRot(:)); minGrndScl = min(grndScale(:));  
maxGrndRot = max(grndRot(:)); maxGrndScl = max(grndScale(:));
stdGrndRot = std(grndRot(:)); stdGrndScl = std(grndScale(:));
fprintf(sprintf('----- Rotation: [%f %f] ----- \n',minGrndRot,maxGrndRot));
fprintf(sprintf('--------- Size: [%f %f] ----- \n',minGrndScl,maxGrndScl));
fprintf(sprintf('----- Rotation Std: %f ----- \n',stdGrndRot));
fprintf(sprintf('--------- Size Std: %f ----- \n',stdGrndScl));

end