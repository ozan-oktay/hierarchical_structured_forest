function generateBarPlot()
% Author: Ozan Oktay
% Date  : 23/10/2015 

close all;
addpath(genpath('/homes/oo2113/workspace/matlab/plotTools/raacampbell13-notBoxPlot-409489c'));
patientId1='14EA01979_ed0_3D'; % might be good for both scale rotation 
patientId2='8SXCUSBX__ed0_3D'; % good for rotation
testing_dataset1 = 'mritestingdata';
testing_dataset2 = 'mribiobankdata';
training_dataset = 'mritrainingdata_sec_large';

% Display the main orientation and scale of the training-data 
rf = 3.14/180;
figure('units','normalized','position',[.5 .1 .2 .5])
trainingdata_pars = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',training_dataset,'/dofs/patientParam.mat'); load(trainingdata_pars);
figure(1);subplot(211), HA=notBoxPlot(z_rot(:)*rf,1,0.6); title('Rotation Distribution (Training Set)','FontSize',14,'fontweight','normal'); d=[HA.data]; set(d(:),'MarkerSize',0.1);
ylabel('Radians','FontSize',14);grid on; hold on; ylim([-50,99]*rf); xlim([0.3,2.4]);
figure(1);subplot(212), HB=notBoxPlot(scale(:),1,0.6); title('Volume Distribution (Training Set)','FontSize',14,'fontweight','normal');    d=[HB.data]; set(d(:),'MarkerSize',0.1);
ylabel('Volume Size Factor','FontSize',14);grid on; hold on; ylim([0.9,1.32]); xlim([0.3,2.4]);


testdata_pars   = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset1,'/dofs/patientParam.mat'); load(testdata_pars);
groun_lm_name_1 = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset1,'/landmarks/',patientId1,'_lm.vtk');
filename        = cellstr(filename); res1=strfind(filename,groun_lm_name_1); nn1=find(~cellfun(@isempty,res1));
rotval1=z_rot(nn1)*rf; scaleval1=scale(nn1);

testdata_pars   = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset2,'/dofs/patientParam.mat'); load(testdata_pars);
groun_lm_name_2 = strcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/',testing_dataset2,'/landmarks/',patientId2,'_lm.vtk');
filename        = cellstr(filename); res2=strfind(filename,groun_lm_name_2); nn2=find(~cellfun(@isempty,res2));
rotval2=z_rot(nn2)*rf; scaleval2=scale(nn2);

figure(1);subplot(211),H1=notBoxPlot(rotval1,1,0.001); d=[H1.data]; set(d(:),'markerfacecolor',[1.0,0.64,0],'color',[0,0.0,0],'MarkerSize',14.0); set(gca,'XTickLabel',[])
figure(1);subplot(211),H2=notBoxPlot(rotval2,1,0.001); d=[H2.data]; set(d(:),'markerfacecolor',[0.4,1,0.4],'color',[0,0.0,0],'MarkerSize',14.0); set(gca,'XTickLabel',[])
h_legend=legend([HA.mu,HA.sdPtch,H1.data,H2.data],'Mean','Std','Image-1','Image-2'); set(h_legend,'FontSize',12); 

figure(1);subplot(212),H1=notBoxPlot(scaleval1,1,0.001);d=[H1.data]; set(d(:),'markerfacecolor',[1.0,0.64,0],'color',[0,0.0,0],'MarkerSize',14.0); set(gca,'XTickLabel',[])
figure(1);subplot(212),H2=notBoxPlot(scaleval2,1,0.001);d=[H2.data]; set(d(:),'markerfacecolor',[0.4,1,0.4],'color',[0,0.0,0],'MarkerSize',14.0); set(gca,'XTickLabel',[])
h_legend=legend([HB.mu,HB.sdPtch,H1.data,H2.data],'Mean','Std','Image-1','Image-2'); set(h_legend,'FontSize',12);



end

