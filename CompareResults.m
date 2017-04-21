clear all; close all;
addpath(genpath(pwd));

subject = 7;
pre_AND_post = true;

if ~ pre_AND_post
    error('only one scan available: cannot compare change');
else
   sheet = 3;
if ~pre_AND_post
   letterOther = letter;
end
if subject < 10
    results_fn_pre  = ['subject0',num2str(subject),'_i_results.mat'];
    results_fn_post = ['subject0',num2str(subject),'_f_results.mat'];
else
    results_fn_pre  = ['subject',num2str(subject),'_i_results.mat'];
    results_fn_post = ['subject',num2str(subject),'_f_results.mat'];
end

% load post scan for comparison
load(results_fn_post)
NslicesPost = Nslices;
sessionOffsetPost = sessionOffset;
upperTopPost =  upperTop;
upperBottomPost =  upperBottom;
upperSlicesPost = upperSlices;
lowerTopPost =  lowerTop;
lowerBottomPost =  lowerBottom;
scanOffsetPost = scanOffset;
diaphramPost = diaphram;
lastSlicePost = lastSlice;
FHpost = FH;

% load pre data 
load(results_fn_pre)
NslicesPre = Nslices;
sessionOffsetPre = sessionOffset;
upperTopPre = upperTop;
upperBottomPre =  upperBottom;
upperSlicesPre = upperSlices;
lowerTopPre =  lowerTop;
lowerBottomPre =  lowerBottom;
scanOffsetPre = scanOffset;
diaphramPre = diaphram;
lastSlicePre = lastSlice;
FHpre = FH; 

first_lower_i_pre = length(upperSlicesPre) + lowerTopPre - upperBottomPre - scanOffsetPre;

overlap_lower_i_pre = (length(upperSlicesPre)+1):(first_lower_i_pre-1);
overlap_upper_i_pre = (length(upperSlicesPre)-length(overlap_lower_i_pre)+1):length(upperSlicesPre);
no_overlap_i_pre    = [1:(length(upperSlicesPre)-length(overlap_lower_i_pre)) first_lower_i_pre:Nslices];

NoverlapPre = length(overlap_lower_i_pre)
% take a weighted average of volumes from both scans for overlaping region
for j = 1:NoverlapPre
    average_weight = j/(NoverlapPre+1);
    SCATvolume(overlap_upper_i_pre(j))   = average_weight*SCATvolume(overlap_lower_i_pre(j))   + (1-average_weight)*SCATvolume(overlap_upper_i_pre(j));
    VATvolume(overlap_upper_i_pre(j))    = average_weight*VATvolume(overlap_lower_i_pre(j))    + (1-average_weight)*VATvolume(overlap_upper_i_pre(j));
    IMATvolume(overlap_upper_i_pre(j))   = average_weight*IMATvolume(overlap_lower_i_pre(j))   + (1-average_weight)*IMATvolume(overlap_upper_i_pre(j));
    organsVolume(overlap_upper_i_pre(j)) = average_weight*organsVolume(overlap_lower_i_pre(j)) + (1-average_weight)*organsVolume(overlap_upper_i_pre(j));
    voidsVolume(overlap_upper_i_pre(j))  = average_weight*voidsVolume(overlap_lower_i_pre(j))  + (1-average_weight)*voidsVolume(overlap_upper_i_pre(j));       
end

% compute shared thorax slices
Nthorax = min(diaphram,diaphramPost)
if diaphramPre > diaphramPost
    thorax_i_pre = (diaphramPre-Nthorax+1):diaphramPre;
else
    thorax_i_pre = 1:diaphramPre;
end

% compute shared abdomen slices (with clipping allowed)
FH_i_pre = length(upperSlices) + lowerTop - FH + 1;
NabdomenPre = FH_i_pre - diaphramPre - NoverlapPre;

FH_i_post = length(upperSlicesPost) + lowerTopPost - FHpost + 1;
NoverlapPost = lowerTopPost - upperBottomPost - scanOffsetPost - 1;
NabdomenPost = FH_i_post - diaphramPost - NoverlapPost;

Nabdomen = min(NabdomenPre,NabdomenPost)
if NabdomenPre > NabdomenPost
    abdomen_i_pre = [(diaphramPre+2):length(upperSlicesPre) first_lower_i_pre:(FH_i_pre-(NabdomenPre-NabdomenPost-1))];
else
    abdomen_i_pre = [(diaphramPre+1):length(upperSlicesPre) first_lower_i_pre:FH_i_pre];
end
length(abdomen_i_pre)

show_i_pre = [thorax_i_pre abdomen_i_pre];
Nshow_pre = length(show_i_pre)

% handle clipped slices whose areas fall outside of the FOV
if lastSlice>0
    Nclipped_post = max((lastSlice-FH),(lastSlicePost-FHpost))
    abdomen_no_clip_i_post =  abdomen_i_pre(1:(end-Nclipped_post));
    abdomen_clipped_i_post =  abdomen_i_pre((end-Nclipped_post+1):end);
end
NabdomenNC_pre = length(abdomen_no_clip_i_post)

% show results
figure, 
imshow(Isaggital(show_i_pre,:),[])
daspect([3 1.40625 1])

% estimate SCAT for clipped slices
for k = 1:Nclipped_post
    SCATbinary  = (SCAT3d(:,:,abdomen_clipped_i_post(k))>0);
    SCATsum1    = sum(SCATbinary);
    SCATsum2    = sum(fliplr(SCATbinary));
    SCATcorrect = sum((SCATsum1==0).*SCATsum2);
    SCATvolume(abdomen_clipped_i_post(k)) = SCATvolume(abdomen_clipped_i_post(k)) + SCATcorrect.*voxelVolumes(abdomen_clipped_i_post(k))/1000.0;
end
aortaVolumePre = aortaVolume;  
PAATvolumePre = PAATvolume; 
heartVolumePre = heartVolume; 
EATvolumePre = EATvolume; 
CATvolumePre = CATvolume; 
voidsVolumePre = voidsVolume; 
lungVolumePre = lungVolume; 
VATvolumePre = VATvolume; 
organsVolumePre = organsVolume; 
IMATvolumePre = IMATvolume; 
SCATvolumePre = SCATvolume;

% load post data
load(results_fn_post)

first_lower_i_post = length(upperSlicesPost) + lowerTopPost - upperBottomPost - scanOffsetPost;

overlap_lower_i_post = (length(upperSlicesPost)+1):(first_lower_i_post-1);
overlap_upper_i_post = (length(upperSlicesPost)-length(overlap_lower_i_post)+1):length(upperSlicesPost);
no_overlap_i_post = [1:(length(upperSlicesPost)-length(overlap_lower_i_post)) first_lower_i_post:NslicesPost];

NoverlapPost = length(overlap_lower_i_post)
% take a weighted average of volumes from both scans for overlaping region
for j = 1:NoverlapPost
    average_weight = j/(NoverlapPost+1);
    SCATvolume(overlap_upper_i_post(j))   = average_weight*SCATvolume(overlap_lower_i_post(j))   + (1-average_weight)*SCATvolume(overlap_upper_i_post(j));
    VATvolume(overlap_upper_i_post(j))    = average_weight*VATvolume(overlap_lower_i_post(j))    + (1-average_weight)*VATvolume(overlap_upper_i_post(j));
    IMATvolume(overlap_upper_i_post(j))   = average_weight*IMATvolume(overlap_lower_i_post(j))   + (1-average_weight)*IMATvolume(overlap_upper_i_post(j));
    organsVolume(overlap_upper_i_post(j)) = average_weight*organsVolume(overlap_lower_i_post(j)) + (1-average_weight)*organsVolume(overlap_upper_i_post(j));
    voidsVolume(overlap_upper_i_post(j))  = average_weight*voidsVolume(overlap_lower_i_post(j))  + (1-average_weight)*voidsVolume(overlap_upper_i_post(j));       
end

% compute shared thorax slices
if diaphramPost > diaphramPre
    thorax_i_post = [(diaphramPost-Nthorax+1):diaphramPost];
else
    thorax_i_post = [1:diaphramPost];
end
length(thorax_i_post)

% compute shared abdomen slices (with clipping allowed)

if NabdomenPost > NabdomenPre
    abdomen_i_post = [(diaphramPost+2):length(upperSlicesPost) first_lower_i_post:(FH_i_post-(NabdomenPost-NabdomenPre-1))];
else
    abdomen_i_post = [(diaphramPost+1):length(upperSlicesPost) first_lower_i_post:FH_i_post];
end

show_i_post = [thorax_i_post abdomen_i_post];
Nshow = length(show_i_post)

% handle clipped slices whose areas fall outside of the FOV
if lastSlice>0
    Nclipped_post = max((lastSlice-FH),(lastSlicePre-FHpre))
    abdomen_no_clip_i_post =  abdomen_i_post(1:(end-Nclipped_post));
    abdomen_clipped_i_post =  abdomen_i_post((end-Nclipped_post+1):end);
end
NabdomenNC_post = length(abdomen_no_clip_i_post)

% show results
figure, 
imshow(Isaggital(show_i_post,:),[])
daspect([3 1.40625 1])

% estimate SCAT for clipped slices
for k = 1:Nclipped_post
    SCATbinary  = (SCAT3d(:,:,abdomen_clipped_i_post(k))>0);
    SCATsum1    = sum(SCATbinary);
    SCATsum2    = sum(fliplr(SCATbinary));
    SCATcorrect = sum((SCATsum1==0).*SCATsum2);
    SCATvolume(abdomen_clipped_i_post(k)) = SCATvolume(abdomen_clipped_i_post(k)) + SCATcorrect.*voxelVolumes(abdomen_clipped_i_post(k))/1000.0;
end

EATvolumePre((Nthorax+1):length(CATvolumePre)) = 0;
EATvolume((Nthorax+1):length(CATvolume)) = 0;
yT = [aortaVolume(show_i_post) - aortaVolumePre(show_i_pre); PAATvolume(show_i_post) - PAATvolumePre(show_i_pre); heartVolume(show_i_post) - heartVolumePre(show_i_pre); EATvolume(show_i_post) - EATvolumePre(show_i_pre); CATvolume(show_i_post) - CATvolumePre(show_i_pre); voidsVolume(show_i_post) - voidsVolumePre(show_i_pre); lungVolume(show_i_post) - lungVolumePre(show_i_pre); VATvolume(show_i_post) - VATvolumePre(show_i_pre); organsVolume(show_i_post) - organsVolumePre(show_i_pre); IMATvolume(show_i_post) - IMATvolumePre(show_i_pre); SCATvolume(show_i_post) - SCATvolumePre(show_i_pre)]'; 
figure, 
bhT = barh(yT,'stacked','LineStyle','none');
set(bhT(1),'FaceColor','magenta');
set(bhT(2),'FaceColor','cyan');
set(bhT(3),'FaceColor','red');
set(bhT(4),'FaceColor',[0.75,0.75,0.75]);
set(bhT(5),'FaceColor',[0,0.75,0.5]);
set(bhT(6),'FaceColor',[0,0,0.5]);
set(bhT(7),'FaceColor','blue');
set(bhT(8),'FaceColor','green');
set(bhT(9),'FaceColor',[0.5,0,0]);
set(bhT(10),'FaceColor',[1 0.667 0.333]);
set(bhT(11),'FaceColor','yellow');
set(gca,'Color','k')
axT = gca; 
axT.YDir = 'reverse';
axT.YLim = [1 Nshow];
xlabel('Volume [mL]')
ylabel('slice number')
title('Total tissue segmentation')

yFp = [PAATvolume(show_i_post) - PAATvolumePre(show_i_pre); EATvolume(show_i_post) - EATvolumePre(show_i_pre); CATvolume(show_i_post) - CATvolumePre(show_i_pre); VATvolume(show_i_post) - VATvolumePre(show_i_pre) + IMATvolume(show_i_post) - IMATvolumePre(show_i_pre); SCATvolume(show_i_post) - SCATvolumePre(show_i_pre)]'; 
figure, 
bhFp = barh(yFp,'stacked','LineStyle','none');
set(bhFp(1),'FaceColor','cyan');
set(bhFp(2),'FaceColor',[0.75,0.75,0.75]);
set(bhFp(3),'FaceColor',[0,0.75,0.5]);
set(bhFp(4),'FaceColor','green');
set(bhFp(5),'FaceColor','yellow');
set(gca,'Color','k')
axFp = gca; 
axFp.YDir = 'reverse';
axFp.YLim = [1 Nshow];
xlabel('Volume [mL]')
ylabel('slice number')
title('Adipose tissue segmentation')

EATvolume(Nthorax+1:length(CATvolume)) = 0;
yF = [PAATvolume(show_i_post) - PAATvolumePre(show_i_pre); EATvolume(show_i_post) - EATvolumePre(show_i_pre); CATvolume(show_i_post) - CATvolumePre(show_i_pre); VATvolume(show_i_post) - VATvolumePre(show_i_pre); IMATvolume(show_i_post) - IMATvolumePre(show_i_pre); SCATvolume(show_i_post) - SCATvolumePre(show_i_pre)]';
figure, 
bhF = barh(yF,'stacked','LineStyle','none');
set(bhF(1),'FaceColor','cyan');
set(bhF(2),'FaceColor',[0.75,0.75,0.75]);
set(bhF(3),'FaceColor',[0,0.75,0.5]);
set(bhF(4),'FaceColor','green');
set(bhF(5),'FaceColor',[1 0.667 0.333]);
set(bhF(6),'FaceColor','yellow');
set(gca,'Color','k')
axF = gca; 
axF.YDir = 'reverse';
axF.YLim = [1 Nshow];
xlabel('Volume [mL]')
ylabel('slice number')
title('Adipose tissue segmentation')

yFreport = [EATvolume(show_i_post) - EATvolumePre(show_i_pre); VATvolume(show_i_post) - VATvolumePre(show_i_pre); SCATvolume(show_i_post) - SCATvolumePre(show_i_pre)]';
figure, 
bhFreport = barh(yFreport,'stacked','LineStyle','none');
set(bhFreport(1),'FaceColor','black');
set(bhFreport(2),'FaceColor','green');
set(bhFreport(3),'FaceColor','yellow');
set(gca,'Color','w')
axFreport = gca; 
axFreport.YDir = 'reverse';
axFreport.YLim = [1 Nshow];
xlabel('Volume [mL]')
ylabel('slice number')
title('Adipose tissue segmentation')

yW = [aortaVolume(show_i_post) - aortaVolumePre(show_i_pre); heartVolume(show_i_post) - heartVolumePre(show_i_pre); organsVolume(show_i_post) - organsVolumePre(show_i_pre)]'; 
figure,
bhW = barh(yW,'stacked','LineStyle','none');
set(bhW(1),'FaceColor','magenta');
set(bhW(2),'FaceColor','red');
set(bhW(3),'FaceColor',[0.5,0,0]);
set(gca,'Color','k')
axW = gca; 
axW.YDir = 'reverse';
axW.YLim = [1 Nshow];
xlabel('Volume [mL]')
ylabel('slice number')
title('Cardiac/Organ segmentation')


yV = [lungVolume(show_i_post) - lungVolumePre(show_i_pre); voidsVolume(show_i_post) - voidsVolumePre(show_i_pre)]'; 
figure, 
bhV = barh(yV,'stacked','LineStyle','none');
set(bhV(1),'FaceColor','blue');
set(bhV(2),'FaceColor',[0,0,0.5]);
set(gca,'Color','k')
axV = gca;
axV.YDir = 'reverse';
axV.YLim = [1 Nshow];
xlabel('Volume [mL]')
ylabel('slice number')
title('Lung/Void segmentation')


yCF = [PAATvolume(show_i_post) - PAATvolumePre(show_i_pre); EATvolume(show_i_post) - EATvolumePre(show_i_pre); CATvolume(show_i_post) - CATvolumePre(show_i_pre)]';
figure, 
bhCF = barh(yCF,'stacked','LineStyle','none');
bhCF(1).FaceColor = 'cyan';
bhCF(2).FaceColor = [0.75,0.75,0.75];
bhCF(3).FaceColor = [0.0,0.75,0.5];
set(gca,'Color','k')
axCF = gca;
axCF.YDir = 'reverse';
axCF.YLim = [1 Nthorax];
xlabel('Volume [mL]')
ylabel('slice number')
title('Cardiac fat segmentation')

yCV = [aortaVolume(show_i_post) - aortaVolumePre(show_i_pre); PAATvolume(show_i_post) - PAATvolumePre(show_i_pre); heartVolume(show_i_post) - heartVolumePre(show_i_pre); EATvolume(show_i_post) - EATvolumePre(show_i_pre); CATvolume(show_i_post) - CATvolumePre(show_i_pre)]'; 
figure, 
bhCV = barh(yCV,'stacked','LineStyle','none');
set(bhCV(1),'FaceColor','magenta');
set(bhCV(2),'FaceColor','cyan');
set(bhCV(3),'FaceColor','red');
set(bhCV(4),'FaceColor',[0.75,0.75,0.75]);
set(bhCV(5),'FaceColor',[0.0,0.75,0.5]);
set(gca,'Color','k')
axCV = gca;
axCV.YDir = 'reverse';
axCV.YLim = [1 Nthorax];
xlabel('Volume [mL]')
ylabel('slice number')
title('Cardiovascular plus fat segmentation')

end