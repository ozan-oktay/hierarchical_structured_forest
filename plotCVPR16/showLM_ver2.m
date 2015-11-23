function showLM_ver2()
% Author: Ozan Oktay
% Date  : 16/11/2015 
close all;
clc;

addpath(genpath('/homes/oo2113/workspace/matlab/plotTools/raacampbell13-notBoxPlot-409489c'));
addpath(genpath('/vol/biomedic/users/oo2113/str_hier_forest_mri'));
testing_dataset    = 'mristacomdata';
training_dataset   = 'mritrainingdata_sec_large';
model              = 'new_prn016_shape_large';
save_dir           = '/homes/oo2113/tmp';

niifilenames       = dir(strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/images/*.nii.gz'));
for loopId=1:numel(niifilenames), niifilename=niifilenames(loopId).name;

    patientId          = strsplit(niifilename,'.nii.gz'); patientId = patientId{1};
    test_image_name    = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/images/',patientId,'.nii.gz');
    strat_lm_name      = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/pems_',model,'/',patientId,'_lm.vtk');
    hough_lm_name      = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/pems/',patientId,'_lm.vtk');

    groun_lm_name      = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/landmarks/',patientId,'.vtk');
    groun_lm_name2     = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/landmarks/',patientId,'_lm.vtk');

    trainingdata_pars  = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',training_dataset,'/dofs/patientParam.mat'); 
    testdata_pars      = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/dofs/patientParam.mat'); 
    poseEstimationsTxt = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/results/poseEstimations_',model,'.txt');

    % Load the nifti image 
    target_image = load_untouch_nii(test_image_name); target_image = target_image.img; 

    % Load the vtk landmark files 
    w2iMat   = world2ImageMat(test_image_name);
    strat_lm = w2iMat * cat(2,vtk2Mat(strat_lm_name),ones(6,1))' +1;
    hough_lm = w2iMat * cat(2,vtk2Mat(hough_lm_name),ones(6,1))' +1;
    groun_lm = w2iMat * cat(2,vtk2Mat(groun_lm_name),ones(6,1))' +1;

    % Show the image and overlay the detected landmarks
    figure('Name','Landmark Detection Results','NumberTitle','off','units','normalized','position',[.1 .1 .4 .6])
    slice = target_image(:,:,round(groun_lm(3,1)));
    imagesc(slice); hold on;
    for i=1:4
        plot(strat_lm(2,i),strat_lm(1,i),'c*','MarkerSize',14);
        plot(hough_lm(2,i),hough_lm(1,i),'r*','MarkerSize',14);
        plot(groun_lm(2,i),groun_lm(1,i),'g*','MarkerSize',14);
    end
    camroll(90)
    colormap('gray');
    h_legend = legend('Stratification Forest','Hough Forest','Manual Annotations','Location','NorthEast');
    set(h_legend,'FontSize',16);
    set(h_legend,'color','k');
    set(h_legend,'TextColor','w')
    set(gca,'YTickLabel',[])
    set(gca,'XTickLabel',[])
    set(gca,'xdir','reverse')

    % Training pose distribution
    load(trainingdata_pars);
    traindata_mean_rot   = mean(z_rot(:));
    traindata_mean_scale = mean(scale(:));

    % Test data ground truth pose information
    load(testdata_pars);
    filename          = cellstr(filename); res=strfind(filename,groun_lm_name2); testPatId=find(~cellfun(@isempty,res));
    testdata_rotval   = z_rot(testPatId); 
    testdata_scaleval = scale(testPatId);

    % Computed pose information 
    firstFileID = fopen(poseEstimationsTxt, 'r');
    formatSpec  = '%s %s %s %f %f %s %f %f %s %f %f';
    firstData   = textscan(firstFileID, formatSpec,'Delimiter',' ');
    for firstPointInd=1:numel(firstData{2})
        fistPointLongName=firstData{2}(firstPointInd);
        firstPointName=strsplit(fistPointLongName{1},'/');
        firstPointName=remSpace(firstPointName{end});
        if (strcmp(firstPointName,strcat(patientId,'.nii.gz')))        
            comp_rot   = firstData{7}(firstPointInd);
            comp_scale = firstData{8}(firstPointInd);
            break;
        end
    end
    
    [r,c,~] = size(target_image);
    t = text(r-20,c-3,sprintf('Subject Id: %s',patientId)); t.FontSize = 12; t.Color='white';

    t = text(r-30,c-3,sprintf('Training Mean Rotation: %2.2f',traindata_mean_rot)); t.FontSize = 12; t.Color='white';
    t = text(r-40,c-3,sprintf('Testing. Data Rotation: %2.2f',testdata_rotval)); t.FontSize = 12; t.Color='white';
    t = text(r-50,c-3,sprintf('Comput Data Rotation: %2.2f',comp_rot)); t.FontSize = 12; t.Color='white';

    t = text(r-70,c-3,sprintf('Training Mean Scale: %2.2f',traindata_mean_scale)); t.FontSize = 12; t.Color='white';
    t = text(r-80,c-3,sprintf('Testing. Data Scale: %2.2f',testdata_scaleval)); t.FontSize = 12; t.Color='white';
    t = text(r-90,c-3,sprintf('Comput Data Scale: %2.2f',comp_scale)); t.FontSize = 12; t.Color='white';
    
    
    figuresavename = strcat(save_dir,'/',patientId,'.png');
    saveas(gcf,figuresavename)
    hold off;
    close all;
   
end
end
function outstr = remSpace(str)
outstr= str(~isspace(str));
end

