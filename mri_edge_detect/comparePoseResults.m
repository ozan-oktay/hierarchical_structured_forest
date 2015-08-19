function comparePoseResults()
%COMPARERESULTS Summary of this function goes here
%   Detailed explanation goes here

firstTxtFile    = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritestingdata/results/poseEstimations_mriSecond_hiera.txt';
workdirectory   = '/vol/biomedic/users/oo2113/str_hier_forest_mri'; addpath(workdirectory);
addpath(genpath(horzcat(workdirectory,'/toolbox')));
  
firstFileID  = fopen(firstTxtFile, 'r');

formatSpec= ['%s %s %s %f %s %f'];
firstData = textscan(firstFileID, formatSpec,'Delimiter',' ');
valueInd  = [4,6];
values    = cell(numel(firstData)/2-1,1);

for firstPointInd=1:numel(firstData{1})
  for lmInd=1:size(values,1), iter=valueInd(lmInd);
    values{lmInd} = [values{lmInd}; firstData{iter}(firstPointInd)];
  end
end

% mean error for pose 
meanErr = sqrt(sum((values{1} - values{2}).^2)/numel(values{1}));

fprintf ('*/*/*/*/*/ Final Results /*/*/*/*/*/ \n');
fprintf (sprintf('\n*/*/*/* total number of images: %d */*/*/*\n',size(values{1},1)));
fprintf ('distances txt name: %s\n',firstTxtFile);
fprintf (sprintf('*/*/*/ Mean Pose Err: %f \n',meanErr));
fprintf ('\n');
fclose(firstFileID);

end
