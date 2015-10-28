function showLM()
% Author: Ozan Oktay
% Date  : 23/10/2015 
close all;
addpath(genpath('/homes/oo2113/workspace/matlab/plotTools/raacampbell13-notBoxPlot-409489c'));
%patientId='8SXCUSBX__ed0_3D'; % might be good for both scale rotation 
patientId='JZNXYS5T__ed0_3D'; % good for rotation
testing_dataset  = 'mribiobankdata';
training_dataset = 'mritrainingdata_sec_large';


test_image_name = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/images/',patientId,'.nii.gz');
strat_lm_name   = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/pems_prn016/',patientId,'_lm.vtk');
hough_lm_name   = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/pems_prn000/',patientId,'_lm.vtk');
groun_lm_name   = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/landmarks/',patientId,'.vtk');

% Load the nifti image 
target_image = load_untouch_nii(test_image_name); target_image = target_image.img; 

% Load the vtk landmark files 
w2iMat   = world2ImageMat(test_image_name);
strat_lm = w2iMat * cat(2,vtk2Mat(strat_lm_name),ones(6,1))' +1;
hough_lm = w2iMat * cat(2,vtk2Mat(hough_lm_name),ones(6,1))' +1;
groun_lm = w2iMat * cat(2,vtk2Mat(groun_lm_name),ones(6,1))' +1;

% Show the image and overlay the detected landmarks
figure('units','normalized','position',[.1 .1 .4 .6])
slice = target_image(:,:,round(groun_lm(3,3)));

    imagesc(slice); hold on;
    for i=1:4
        plot(strat_lm(2,i),strat_lm(1,i),'c*','MarkerSize',16);
        plot(hough_lm(2,i),hough_lm(1,i),'r*','MarkerSize',16);
        plot(groun_lm(2,i),groun_lm(1,i),'g*','MarkerSize',16);
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

hold off;

% Display the main orientation and scale of the training-data 
figure('units','normalized','position',[.5 .1 .3 .4])
trainingdata_pars = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',training_dataset,'/dofs/patientParam.mat'); load(trainingdata_pars);
figure(2);subplot(121), H=notBoxPlot(z_rot(:),1,0.2); title('Rotation Dist of the TrainingSet','FontSize',10); d=[H.data]; set(d(:),'MarkerSize',0.1);
ylabel('Degrees wrt Reference Image','FontSize',12);grid on; hold on; ylim([-50,99]);
figure(2);subplot(122), H=notBoxPlot(scale(:),1,0.2); title('Volume Dist of the TrainingSet','FontSize',10);    d=[H.data]; set(d(:),'MarkerSize',0.1);
ylabel('Volume Scale Factor','FontSize',12);grid on; hold on; ylim([0.9,1.32])

testdata_pars = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/dofs/patientParam.mat'); load(testdata_pars);
groun_lm_name = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset,'/landmarks/',patientId,'_lm.vtk');
filename      = cellstr(filename); res=strfind(filename,groun_lm_name); nn=find(~cellfun(@isempty,res));
rotval=z_rot(nn); scaleval=scale(nn);
figure(2);subplot(121),H=notBoxPlot(rotval,1,0.001); d=[H.data]; set(d(:),'markerfacecolor',[0.4,1,0.4],'color',[0,0.4,0],'MarkerSize',15.0); set(gca,'XTickLabel',[])
h_legend=legend([H.data],'Testing Image'); set(h_legend,'FontSize',10);
figure(2);subplot(122),H=notBoxPlot(scaleval,1,0.001);d=[H.data]; set(d(:),'markerfacecolor',[0.4,1,0.4],'color',[0,0.4,0],'MarkerSize',15.0); set(gca,'XTickLabel',[])
h_legend=legend([H.data],'Testing Image'); set(h_legend,'FontSize',10);

end

