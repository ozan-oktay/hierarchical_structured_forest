function biobank_correct_orientation_pems
%CORRECTORIENTATIONOFBIOPEMS Summary of this function goes here
%   Detailed explanation goes here

pems_directory    = '/vol/biomedic/users/oo2113/str_hier_forest_mri/biobankdata/pems_reg';
biobank_directory = '/vol/vipdata/data/biobank/cardiac/nifti';
output_directory  = '/vol/biomedic/users/oo2113/for_Wenjia';

% Find the generated pem files 
pems_names = dir(strcat(pems_directory,'/*.nii.gz')); 
pems_names = {pems_names(:).name}'; 
pems_names(ismember(pems_names,{'.','..'})) = [];

% Find biobank folders 
d = dir(strcat(biobank_directory,'/*')); 
biobankfoldernames = {d([d(:).isdir]).name}'; 
biobankfoldernames(ismember(biobankfoldernames,{'.','..'})) = [];

for pems_Ind=1:numel(pems_names)
  
  pems_name    = pems_names{pems_Ind};
  patient_name = strsplit(pems_name,'_'); 
  patient_name = horzcat(patient_name{1});
  
  for bio_Ind=1:numel(biobankfoldernames)
    
    biofolder_name = biobankfoldernames{bio_Ind};
    biofolder_name = strsplit(biofolder_name,'_'); 
    biofolder_name = horzcat(biofolder_name{1});
    
    if isequal(patient_name,biofolder_name)
      
      targetname   = horzcat(biobank_directory,'/',biobankfoldernames{bio_Ind},'/cine_sax.nii.gz');
      inputlmname  = horzcat(pems_directory   ,'/',pems_names{pems_Ind});
      outputlmname = horzcat(output_directory ,'/',pems_names{pems_Ind});
      tmpimagename = strcat(tempname('/homes/oo2113/tmp'),'.nii.gz');
      
      system(horzcat('headertool',' ',inputlmname,' ',tmpimagename,' ','-reset'));
      system(horzcat('headertool',' ',tmpimagename,' ',outputlmname,' ','-orientation',' ','1 0 0 0 -1 0 0 0 -1'));
      system(horzcat('transformation',' ',outputlmname,' ',outputlmname,' ','-target',' ',tmpimagename,' ','-linear -matchInputType'));
      system(horzcat('headertool',' ',outputlmname,' ',outputlmname,' ','-targetOriginAndOrient',' ',targetname));
      
      delete(tmpimagename);
    end
  end
end
end

