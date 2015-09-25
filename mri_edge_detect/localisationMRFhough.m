%% random forest hough appearance model for finding landmarks
clear all; close all;
testnum=1; %select any test case
hough=load_nii(['~/Downloads/mristacomdata/pems/cSAX_',num2str(testnum),'_ed_lm.nii.gz']); hough=hough.img;
%% generate a mean location model for all landmarks
model_test=load(['/Users/mattias/Downloads/mristacomdata/landmarks/cSAX_',num2str(testnum),'_ed.txt']);
model=zeros(6,3,8); trainnum=[1,2,4:10]; trainnum(testnum)=[]; %leave out the test locations
for i=1:length(trainnum);
    model1=load(['/Users/mattias/Downloads/mristacomdata/landmarks/cSAX_',num2str(trainnum(i)),'_ed.txt']);
    model(:,:,i)=model1;
end
mean_model=mean(model,3);
centre_gt=model_test(:,:)+1;

%% perform MRF regularisation
%both following parameters have to be emperically chosen
alpha=5; % a greater alpha gives weaker regularisation
sigma=.3; % fall-off for exponential function for unary cost

% generate unary term based on hough image
match_dist_all=zeros(size(hough),'single');
hough_n=hough./max(hough(:)); %scale max = 1
for l=1:6;
    match_dist_all(:,:,:,l)=exp(-(hough_n(:,:,:,l))./sigma);
end

% assume simple star-graph model: every landmark influences every other
% landmark, but no joint optimisation is done
for Label=1:6;
    regular_dist=match_dist_all(:,:,:,Label);
    h=hough(:,:,:,Label);
    %error without regularisation (in mm)
    [val,ind]=max(h(:)); [y_m,x_m,z_m]=ind2sub(size(regular_dist),ind);
    locErrBefore(Label)=sqrt((1.25*centre_gt(Label,1)-1.25*y_m).^2+(1.25*centre_gt(Label,2)-1.25*x_m).^2+(2*centre_gt(Label,3)-2*z_m).^2);
    %weigh the edges in the graph based on their inverse squared Euclidean
    %distance (divided by sum/5 for normalisation)
    length_model=(sum((repmat(mean_model(Label,:),6,1)-mean_model).^2,2));
    length_model=length_model./(sum(length_model)/5);
    for l=1:6;
        if(Label~=l)
            offset=single(mean_model(Label,:)-mean_model(l,:));
            %parts-based model with squared regularisation term
            reg_dist_aux=distTrans3doff(match_dist_all(:,:,:,l),single(offset),alpha/length_model(l));
            %simply add the infered pairwise cost-maps to unary term
            regular_dist=regular_dist+reg_dist_aux;
        end
    end
    % eval optimum and compare to ground truth after regularisation
    [val,ind]=min(regular_dist(:));
    [y_r,x_r,z_r]=ind2sub(size(regular_dist),ind);
    locErrReg(Label)=sqrt((1.25*centre_gt(Label,1)-1.25*y_r).^2+(1.25*centre_gt(Label,2)-1.25*x_r).^2+(2*centre_gt(Label,3)-2*z_r).^2);
end
[mean(locErrBefore),mean(locErrReg)]

figure; plot([locErrBefore',locErrReg']);
