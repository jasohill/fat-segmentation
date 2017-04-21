%------------------------------------------------------------------------
% Author: Dr. Jason E. Hill, post-doctoral researcher
% Institution: CNG at TTU in a working association with TTNI
% Date: 9 MAR 2016
% Updated: 16 MAR 2016
%------------------------------------------------------------------------

clear all, close all, clc;

addpath(genpath(pwd))
addpath('../../auxillary')
addpath('../../auxillary/NIfTI_20140122')
addpath('C:/spm/spm8')

% subject: MF0202
% session: 150522
% date:    12 OCT 2015


%% scan info

NIFTI_file_name_F_upper = '20151116_150522t1vibedixontrap4bh320s004a1001.nii.gz';
NIFTI_file_name_W_upper = '20151116_150522t1vibedixontrap4bh320s005a1001.nii.gz';
NIFTI_file_name_F_lower = '20151116_150522t1vibedixontrap4bh320s009a1001.nii.gz';
NIFTI_file_name_W_lower = '20151116_150522t1vibedixontrap4bh320s010a1001.nii.gz';

NIFTI_file_name_F_upper_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_F_0004/20151116_150522t1vibedixontrap4bh320s004a1001.nii';
NIFTI_file_name_W_upper_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_W_0005/20151116_150522t1vibedixontrap4bh320s005a1001.nii';
NIFTI_file_name_F_lower_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_F_0009/20151116_150522t1vibedixontrap4bh320s009a1001.nii';
NIFTI_file_name_W_lower_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_W_0010/20151116_150522t1vibedixontrap4bh320s010a1001.nii';

NIFTI_file_name_F_upper_reslice = '20151116_150522t1vibedixontrap4bh320s004a1001_resliced.nii';
NIFTI_file_name_W_upper_reslice = '20151116_150522t1vibedixontrap4bh320s005a1001_resliced.nii';
NIFTI_file_name_F_lower_reslice = '20151116_150522t1vibedixontrap4bh320s009a1001_resliced.nii';
NIFTI_file_name_W_lower_reslice = '20151116_150522t1vibedixontrap4bh320s010a1001_resliced.nii';

NIFTI_file_name_F_upper_reslice2 = '20151116_150522t1vibedixontrap4bh320s004a1001_resliced_resliced.nii';
NIFTI_file_name_W_upper_reslice2 = '20151116_150522t1vibedixontrap4bh320s005a1001_resliced_resliced.nii';
NIFTI_file_name_F_lower_reslice2 = '20151116_150522t1vibedixontrap4bh320s009a1001_resliced_resliced.nii';
NIFTI_file_name_W_lower_reslice2 = '20151116_150522t1vibedixontrap4bh320s010a1001_resliced_resliced.nii';

voxelSize_up    = [1.40625, 1.40625, 2.5];
interslice_spacing_fraction = 0.2;
voxelVolume_up  = prod(voxelSize_up)*(1+interslice_spacing_fraction);
voxelSize_lo    = [1.40625, 1.40625, 2.5];
voxelVolume_lo  = prod(voxelSize_up)*(1+interslice_spacing_fraction);
% voxelSize_up(3) = (1+interslice_spacing_fraction)*voxelSize_up(3);
% voxelSize_lo(3) = (1+interslice_spacing_fraction)*voxelSize_up(3);

% matching slices for stitching purposes
showSliceUp = 22;
showSliceLo = 112;

%% Upper scan rectification

% load fat dominated signal data for upper scan
niiFup = load_nii(NIFTI_file_name_F_upper);
niiFupMax = max(max(max(niiFup.img)));
% load water selected data for upper scan
niiWup = load_nii(NIFTI_file_name_W_upper);
niiWupMax = max(max(max(niiWup.img)));

% make threshold image from both data sets
Bup = Make_DIXON_threshold_image(niiFup.img,niiWup.img);

BupCoM = center_of_mass(Bup>0.5)
BupShowSliceCoM = center_of_mass(Bup(:,:,showSliceUp)>0.5)

[tiltAngleUp,rotAngleUp] = Find_table_tilt_angle(Bup, voxelSize_up,[21,92],[],false); %[21,92]; [1,65] %true to see table lines
tiltAngleUp
rotAngleUp

addpath(genpath('C:\spm\spm8'))

scanUp.voxel.size = voxelSize_up;
scanUp.voxel.dimensions = size(niiFup.img);
scanUp.voxel.spacing = [0 0 0];    
scanUp.tiltAngle = tiltAngleUp;

V_F_up = STANCE_load_volume(NIFTI_file_name_F_upper);

V_F_up.mat

tic
[V_F_up_reslice,Y_F_up_reslice] = STANCE_reslice_volume(V_F_up,scanUp,[],true,true,3);
toc

V_F_up_reslice_bottom_Y = (2*V_F_up.mat(2,4) - V_F_up_reslice.mat(2,4))

niiFup_reslice = load_nii(NIFTI_file_name_F_upper_reslice);
figure, imshow(fliplr(rot90(niiFup_reslice.img(:,:,showSliceUp))),[]);

% [x_cm_Up, angleL_Up, angleR_Up] = Find_plane_angle(Bup, voxelSize_up, 2.5,[28,60],[],true); %[11,66] %[5,71]
% x_cm_Up
% angleL_Up
% angleR_Up
% % figure, imshow(Bup(:,:,5),[]);
% % figure, imshow(Bup(:,:,30),[]);
% % figure, imshow(Bup(:,:,50),[]);
% % figure, imshow(Bup(:,:,65),[]);
% [x_cm_Ups, y_cm_Ups, angleXZ_Up, angleYZ_Up] = Find_cm_angle(Bup, voxelSize_up, 2.5,[21,92],[],true); %[11,66] %[5,71]
% angleXZ_Up
% %angleYZ_Up

%translationUp = [(160+V_F_up_reslice.mat(1,4)/V_F_up_reslice.mat(1,1)),...
%    ((130.5-mean(x_cm_Up))/V_F_up_reslice.mat(2,2)),0]

% repeat for water saturated data
V_W_up = STANCE_load_volume(NIFTI_file_name_W_upper);
tic
[V_W_up_reslice,Y_W_up_reslice] = STANCE_reslice_volume(V_W_up,scanUp,[],true,true,3);
toc

niiWup_reslice = load_nii(NIFTI_file_name_W_upper_reslice);
%figure, imshow(niiWup_reslice.img(:,:,showSliceUp),[]);
Bup_reslice = Make_DIXON_threshold_image(niiFup_reslice.img,niiWup_reslice.img);

BupCoM = center_of_mass(Bup_reslice>0.5)
BupShowSliceCoM = center_of_mass(Bup_reslice(:,:,showSliceUp)>0.5)

translationUp = [BupShowSliceCoM(1)-261/2.0,...
                 BupShowSliceCoM(2)-321/2.0,0]

%figure, imshow(flipud(Y_F_up_reslice(:,:,1)),[]);
%figure, imshow(flipud(Y_F_up_reslice(:,:,27)),[]);
%figure, imshow(flipud(Y_F_up_reslice(:,:,60)),[]);
%figure, imshow(imrotate(flipud(Y_F_up_reslice(:,:,showSliceUp)),90),[]);
%figure, imshow(flipud(Y_F_up_reslice(:,:,120)),[]);

V_F_up_reslice.mat

tic
Y_F_up_reslice_trans = imtranslate(Y_F_up_reslice,translationUp,'cubic','OutputView','same');
toc

tic
Y_W_up_reslice_trans = imtranslate(Y_W_up_reslice,translationUp,'cubic','OutputView','same');
toc

%figure, imshow(flipud(Y_F_up_reslice_trans(:,:,1)),[]);
%figure, imshow(flipud(Y_F_up_reslice_trans(:,:,27)),[]);
%figure, imshow(flipud(Y_F_up_reslice_trans(:,:,60)),[]);
figure, imshow(fliplr(rot90(flipud(Y_F_up_reslice_trans(:,:,showSliceUp)))),[]);
%figure, imshow(flipud(Y_F_up_reslice_trans(:,:,120)),[]);

tic
for z = 1:120
    Y_F_up_reslice_trans_rot(:,:,z) = imrotate(Y_F_up_reslice_trans(:,:,z),-rotAngleUp,'bicubic','crop');
end
toc
Y_F_up_reslice_trans_rot(Y_F_up_reslice_trans_rot<0) = 0.0;
Y_F_up_reslice_trans_rot(Y_F_up_reslice_trans_rot>niiFupMax) = niiFupMax;

V_F_up_reslice.mat(1,4) = sign(V_F_up_reslice.mat(1,4))*(abs(V_F_up_reslice.mat(1,4)) + translationUp(1)*V_F_up_reslice.mat(1,1));
V_F_up_reslice.mat(2,4) = sign(V_F_up_reslice.mat(2,4))*(abs(V_F_up_reslice.mat(2,4)) + translationUp(2)*V_F_up_reslice.mat(2,2));

V_F_up_reslice = spm_write_vol(V_F_up_reslice,Y_F_up_reslice_trans_rot);

% for z = 1:120
%     Y_F_up_reslice_trans_rot(:,:,z) = imrotate(Y_F_up_reslice_trans(:,:,z),-rotAngleUp,'bilinear','crop')';
% end
% temp_dim = V_F_up_reslice.dim;
% V_F_up_reslice.dim(1) = temp_dim(2);
% V_F_up_reslice.dim(2) = temp_dim(1);
% 
% %figure, imshow(flipud(Y_F_up_reslice_trans_rot(:,:,showSliceUp)),[]);
% 
% V_F_up_reslice.mat(1,4) = sign(V_F_up_reslice.mat(1,4))*(abs(V_F_up_reslice.mat(1,4)) + translationUp(1)*V_F_up_reslice.mat(1,1));
% V_F_up_reslice.mat(2,4) = sign(V_F_up_reslice.mat(2,4))*(abs(V_F_up_reslice.mat(2,4)) + translationUp(2)*V_F_up_reslice.mat(2,2));
% 
% V_F_up_reslice = spm_write_vol(V_F_up_reslice,Y_F_up_reslice_trans_rot);

% plane_angle_Up = 0.25*(angleL_Up + angleR_Up + angleXZ_Up)
% 
% scanUp2 = scanUp;
% scanUp2.voxel.size = [voxelSize_up(2) voxelSize_up(1) voxelSize_up(3)];
% scanUp2.voxel.dimensions = [size(niiFup.img,2) size(niiFup.img,1) size(niiFup.img,3)];   
% scanUp2.tiltAngle = plane_angle_Up;
% 
% [V_F_up_reslice,Y_F_up_reslice] = STANCE_reslice_volume(V_F_up_reslice,scanUp2,[],true,true);
% 
% for z = 1:120
%     Y_F_up_reslice_trans(:,:,z) = Y_F_up_reslice(:,:,z)';
% end
% V_F_up_reslice.dim = temp_dim;
% 
% V_F_up_reslice = spm_write_vol(V_F_up_reslice,Y_F_up_reslice_trans);

Bup_reslice = Make_DIXON_threshold_image(niiFup_reslice.img,niiWup_reslice.img);
[tiltAngleUp_reslice,rotAngleUp_reslice] = Find_table_tilt_angle(Bup_reslice, voxelSize_up,[21,92],[],false); %[21,92] %[20,100] %[7,65]
tiltAngleUp_reslice
rotAngleUp_reslice

Y_W_up_reslice_trans = imtranslate(Y_W_up_reslice,translationUp,'cubic','OutputView','same');
tic
for z = 1:120
    Y_W_up_reslice_trans_rot(:,:,z) = imrotate(Y_W_up_reslice_trans(:,:,z),-rotAngleUp,'bicubic','crop');
end
toc
Y_W_up_reslice_trans_rot(Y_W_up_reslice_trans_rot<0) = 0.0;
Y_W_up_reslice_trans_rot(Y_W_up_reslice_trans_rot>niiWupMax) = niiWupMax;

V_W_up_reslice.mat(1,4) = V_F_up_reslice.mat(1,4);
V_W_up_reslice.mat(2,4) = V_F_up_reslice.mat(2,4);

V_W_up_reslice = spm_write_vol(V_W_up_reslice,Y_W_up_reslice_trans_rot);

% Y_W_up_reslice_trans = imtranslate(Y_W_up_reslice,translationUp,'linear','OutputView','same');
% for z = 1:120
%     Y_W_up_reslice_trans_rot(:,:,z) = imrotate(Y_W_up_reslice_trans(:,:,z),-rotAngleUp,'bilinear','crop')';
% end
% V_W_up_reslice.dim(1) = temp_dim(2);
% V_W_up_reslice.dim(2) = temp_dim(1);
% V_W_up_reslice.mat(1,4) = V_F_up_reslice.mat(1,4);
% V_W_up_reslice.mat(2,4) = V_F_up_reslice.mat(2,4);
% 
% V_W_up_reslice = spm_write_vol(V_W_up_reslice,Y_W_up_reslice_trans_rot);
% 
% % [V_W_up_reslice,Y_W_up_reslice] = STANCE_reslice_volume(V_W_up_reslice,scanUp2,[],true,true);
% % 
% for z = 1:120
%     Y_W_up_reslice_trans(:,:,z) = Y_W_up_reslice(:,:,z)';
% end
% V_W_up_reslice.dim = temp_dim;
% V_W_up_reslice = spm_write_vol(V_W_up_reslice,Y_W_up_reslice_trans);


% Bup_reslice = Make_DIXON_threshold_image(niiFup_reslice.img,niiWup_reslice.img);
% 
% [tiltAngleUp_reslice,rotAngleUp_reslice] = Find_table_tilt_angle(Bup_reslice, voxelSize_up,[21,92],[],false); %[20,100] %[7,65]
% tiltAngleUp_reslice
% rotAngleUp_reslice

% [x_cm_Up, angleL_Up, angleR_Up] = Find_plane_angle(Bup_reslice, voxelSize_up, 2.5,[28,60],[],true);%[5,71]
% x_cm_Up
% angleL_Up
% angleR_Up
% 
% [x_cm_Ups, y_cm_Ups, angleXZ_Up, angleYZ_Up] = Find_cm_angle(Bup_reslice, voxelSize_up, 2.5,[21,92],[],true); %[5,71]
% angleXZ_Up
% %angleYZ_Up
% 
% plane_angle_Up_reslice = 0.25*(angleL_Up + angleR_Up + angleXZ_Up)

V_W_up_reslice.mat

% load rectified data

niiWup_reslice = load_nii(NIFTI_file_name_W_upper_reslice);
figure, imshow(fliplr(rot90(niiWup_reslice.img(:,:,showSliceUp))),[]);

niiFup_reslice = load_nii(NIFTI_file_name_F_upper_reslice);
figure, imshow(fliplr(rot90(niiFup_reslice.img(:,:,showSliceUp))),[]);

%% register landmarks (key anatomical features)
sessionOffset  = 7;

upperBottom = 10;
Z_L3_L4u    = 3;
Z_L2_L3u    = 15;
Z_L1_L2u    = 26;  % middle of L1-L2 intervertebral disk slice (upper scan)
% Z_L1b        = 29;  % bottom of L1 vertebra slice
Z_T4_T5t    = 113; % top of T4-T5 intervertebral disk slice
%upperSlices = Z_T4_T5t:-1:Z_L1_L2u;
Z_T9_T10t      = 59;  % top of T10 vertebra disk slice (use for bottom of thoracic cavity) 
% NOTE: diaphraim is at slices = 60-68
% - register from the bottom of T9 (slice 54)
upperTop = 113;

% NOTE: diaphraim is at slices = 60-68
% - register from the bottom of T9 (slice 54)
heartMax  = upperTop-76+1;
heartApex = upperTop-63+1;
diaphram  = upperTop-Z_T9_T10t+1; % used as bottom of lung/thorax

upperSlices = upperTop:-1:upperBottom;

lowerBottom = 4;
% PF       = 14;   % use the bottom of FH sphere as the pelvic floor
lowerOffset = 15;
FH       = 22;  % widest slice of the formoral head
lastSlice = 53;  % i = 161, 47th 63rd slice from top
Z_L5_S1b = 60;  % bottom of L5-S1 intervertebral disk slice
Z_L4_L5  = 79;   % middle of L4-L5 intervertebral disk slice
Z_L3_L4  = 93;   % middle of L3-L4 intervertebral disk slice
UM       = 96;   % umbilicus (belly button)
Z_L2_L3  = 105;  % middle of L2-L3 intervertebral disk slice
Z_L1_L2  = 116;  % middle of L1-L2 intervertebral disk slice
%Z_L1_L2t    = 118; % top of L2 vertebra disk
lowerTop = 118;
lowerSlices = lowerTop:-1:lowerBottom;
% Z_L1_L2t = 113;  % top of L1-L2 intervertebral disk
%lowerSlices = Z_L1_L2:-1:FH+2; %FH+lowerOffset;

% correct scan offset
scanOffset = 90  %Z_L3_L4 - Z_L3_L4u

first_lower_l = 121 - lowerTop + upperBottom + scanOffset;

% last slice without clipping in lower scan: slice # 56, (i = 159)

Nslices = length(upperSlices) + length(lowerSlices);% + repeatEstimators;
slices = [upperSlices, lowerSlices]; %, Z_L2_L3u-2:Z_L2_L3u+2]; % FH-2:FH+2, Z_L2_L3u-2:Z_L2_L3u+2];
levels(1:length(upperSlices))         = 1;
levels(length(upperSlices)+1:Nslices) = 2;
voxelVolumes(1:length(upperSlices))         = voxelVolume_up;
voxelVolumes(length(upperSlices)+1:Nslices) = voxelVolume_lo;

% NOTE: belly button is from i = 108-114
bellybutton = (length(upperSlices) + lowerTop - UM + 1 - 6):(length(upperSlices) + lowerTop - UM + 1 + 5);
bbx = 70;

% parameter settings
cutLineLength = 82; %40;

backgroundThreshold   = 10; %15%10
vetoFactorF = 4.5;
vetoFactorW = 3.5;
foregroundThreshold1  = 74; %74 %59
foregroundThreshold2  = 10;
foregroundThresholdF2 = 145;
foregroundThresholdV  = 175; %103;
foregroundThresholdV2 = 230;
foregroundThresholdW  = 40;

nbrThresholdsF = 4; %4; %7;
nbrThresholdsW = 7; %4; %9;
nbrThresholdsV = 7; %4; %7;

seSCAT = strel('rectangle',[7,17]);  
VATprelim      = 17;     % opening for preliminary VAT segmentation

dims = size(niiFup.img);
% correction of inner mask by removing effect of mammaries
innerVeto = zeros(1,dims(3));
innerVeto(114:120) = 261 - 181;
innerVeto(113)     = 261 - 182;
innerVeto(112)     = 261 - 183;
innerVeto(111)     = 261 - 184;
innerVeto(109:110) = 261 - 185;
innerVeto(107:108) = 261 - 186;
innerVeto(105:106) = 261 - 187;
innerVeto(103:104) = 261 - 188;
innerVeto(101:102) = 261 - 189;
innerVeto(98:100)  = 261 - 190;
innerVeto(93:97)   = 261 - 191;
innerVeto(88:92)   = 261 - 192;
innerVeto(86:87)   = 261 - 193;
innerVeto(82:85)   = 261 - 194;
innerVeto(77:81)   = 261 - 195;
innerVeto(75:76)   = 261 - 196;
innerVeto(73:74)   = 261 - 197;
innerVeto(70:72)   = 261 - 198;
innerVeto(68:69)   = 261 - 199;
innerVeto(65:67)   = 261 - 200;
%innerVeto      = 70;     % used for correcting the body cavity (inner) mask
%innerVetoLimit = 35;
voidThreshold  = 5000;   % threshold of count of number of void pixels
lungErode      = 6; %6

aortaThreshold = 80;     % threshold of 8-bit water dominated signal data
aortaMinArea   = 70;
aortaDistance  = 50;
aortaRadius    = 9;

notLungDilate  = 45; %45%42%21;
innerDilate    = 12; %13;
heartShift     = 25;
heartDilate    = 45; %45 %42
heartAreaThreshold = 400;
CATareaThreshold = 64;%160;

CATdilate      = 12; %7; %5; %3; %4; %6
CATmargin      = 3; %6;
PAATdilate     = 10; %6;

%CATdilateCorrect = zeros(1,Nslices);
%CATdilateCorrect(heartMax-25) = +2;
%CATdilateCorrect(heartMax-24:heartMax-19) = +3; %+5;
%CATdilateCorrect(heartMax-18) = +2;
%CATdilateCorrect(heartMax) = +2;
%CATdilateCorrect(heartMax+1:heartMax+3) = +3;
%CATdilateCorrect(heartMax+4) = +2;
%CATdilateCorrect(heartApex-4:heartApex) = +3;

% initializations
SCATvolume   = zeros(1,Nslices); % [mL]
VATvolume    = zeros(1,Nslices); % [mL]
organsVolume = zeros(1,Nslices); % [mL]
voidsVolume  = zeros(1,Nslices); % [mL]
lungVolume   = zeros(1,Nslices); % [mL]
heartVolume  = zeros(1,Nslices); % [mL]
aortaVolume  = zeros(1,Nslices); % [mL]
CATvolume    = zeros(1,Nslices); % [mL]
PAATvolume   = zeros(1,Nslices); % [mL]

% define aorta seed in native coordinates of level 1 (upper scan)
aortaSeed = false(dims(2),dims(1),dims(3));
% vertical section
aortaSeed(261-120,165,40)      = 1;
aortaSeed(261-119,165,41)      = 1;
aortaSeed(261-118,166,42)      = 1;
aortaSeed(261-117,166,43)      = 1;
aortaSeed(261-116,166,44:46)   = 1;
aortaSeed(261-115,166,47:59)   = 1;
aortaSeed(261-113,166,60:63)   = 1;
aortaSeed(261-113,165,64:65)   = 1;
aortaSeed(261-112,166,66:69)   = 1;
aortaSeed(261-112,167,70:75)   = 1;
aortaSeed(261-111,167,76:79)   = 1;
aortaSeed(261-110,168,80:85)   = 1;
aortaSeed(261-110,168,85:87)   = 1;
aortaSeed(261-110,169,88:89)   = 1;
aortaSeed(261-110,169,90:91)   = 1;
aortaSeed(261-110,170,92:98)   = 1; 
aortaSeed(261-110,171,99:102)  = 1;
aortaSeed(261-111,171,103:105) = 1;
aortaSeed(261-112,172,106:107) = 1;
aortaSeed(261-113,172,108)     = 1;
aortaSeed(261-113,171,109)     = 1;
aortaSeed(261-113,174,110:112) = 1;
% curved section
aortaSeed(261-144,167,109:110)      = 1;
aortaSeed(261-144,168,111)          = 1;
aortaSeed(261-121,175,109:upperTop) = 1;
aortaSeed(261-127,175,109:upperTop) = 1;
aortaSeed(261-130,174,109:upperTop) = 1;
aortaSeed(261-134,173,109:upperTop) = 1;
aortaSeed(261-150,172,109:upperTop) = 1;
aortaSeed(261-141,169,109:upperTop) = 1;


%% lower scan

% load fat dominated signal data for lower scan
niiFlo = load_nii(NIFTI_file_name_F_lower);
niiFloMax = max(max(max(niiFlo.img)))
% load water selected data for lower scan
niiWlo = load_nii(NIFTI_file_name_W_lower);
niiWloMax = max(max(max(niiWlo.img)))

% make threshold image from both data sets
Blo = Make_DIXON_threshold_image(niiFlo.img,niiWlo.img);

BloCoM = center_of_mass(Blo>0.5)
BloShowSliceCoM = center_of_mass(Blo(:,:,showSliceLo)>0.5)

[tiltAngleLo,rotAngleLo] = Find_table_tilt_angle(Blo, voxelSize_lo,[5,50],[],false); %[20,60] %60 % [20,first_lower_l]
tiltAngleLo
rotAngleLo

% [x_cm_Lo, angleL_Lo, angleR_Lo] = Find_plane_angle(Blo, voxelSize_lo, 2.5,[55,115],[],true);
% x_cm_Lo
% angleL_Lo
% angleR_Lo
% [x_cm_Los, y_cm_Los, angleXZ_Lo, angleYZ_Lo] = Find_cm_angle(Blo, voxelSize_up, 2.5,[55,115],[],true);
% angleXZ_Lo
% %angleYZ_Lo

scanLo.voxel.size = voxelSize_lo;
scanLo.voxel.dimensions = size(niiFlo.img);
scanLo.voxel.spacing = [0 0 0];    
scanLo.tiltAngle = tiltAngleLo;

V_F_lo = STANCE_load_volume(NIFTI_file_name_F_lower);

V_F_lo.mat

[V_F_lo_reslice,Y_F_lo_reslice] = STANCE_reslice_volume(V_F_lo,scanLo,[],true,true);

V_F_lo_reslice.mat

niiFlo_reslice = load_nii(NIFTI_file_name_F_lower_reslice);
figure, imshow(fliplr(rot90(niiFlo_reslice.img(:,:,showSliceLo))),[]);

% repeat for water saturated data
V_W_lo = STANCE_load_volume(NIFTI_file_name_W_lower);
[V_W_lo_reslice,Y_W_lo_reslice] = STANCE_reslice_volume(V_W_lo,scanLo,[],true,true,3);

niiWlo_reslice = load_nii(NIFTI_file_name_W_lower_reslice);
%figure, imshow(niiWlo_reslice.img(:,:,showSliceLo),[]);
Blo_reslice = Make_DIXON_threshold_image(niiFlo_reslice.img,niiWlo_reslice.img);

BloCoM = center_of_mass(Blo_reslice>0.5)
BloShowSliceCoM = center_of_mass(Blo_reslice(:,:,showSliceLo)>0.5)

%translationLo = [(160+V_F_lo_reslice.mat(1,4)/V_F_lo_reslice.mat(1,1)),...
%    (130+V_F_lo_reslice.mat(2,4)/V_F_lo_reslice.mat(2,2)),0]
translationLo = [BloShowSliceCoM(1)-261/2.0,...
                 BloShowSliceCoM(2)-321/2.0,0]

%figure, imshow(flipud(Y_F_lo_reslice(:,:,1)),[]);
%figure, imshow(flipud(Y_F_lo_reslice(:,:,27)),[]);
%figure, imshow(flipud(Y_F_lo_reslice(:,:,60)),[]);
%figure, imshow(flipud(Y_F_lo_reslice(:,:,showSliceLo)),[]);
%figure, imshow(flipud(Y_F_lo_reslice(:,:,120)),[]);

Y_F_lo_reslice_trans = imtranslate(Y_F_lo_reslice,translationLo,'cubic','OutputView','same');

%figure, imshow(flipud(Y_F_lo_reslice_trans(:,:,1)),[]);
%figure, imshow(flipud(Y_F_lo_reslice_trans(:,:,27)),[]);
%figure, imshow(flipud(Y_F_lo_reslice_trans(:,:,60)),[]);
%figure, imshow(flipud(Y_F_lo_reslice_trans(:,:,showSliceLo)),[]);
%figure, imshow(flipud(Y_F_lo_reslice_trans(:,:,120)),[]);
for z = 1:120
    Y_F_lo_reslice_trans_rot(:,:,z) = imrotate(Y_F_lo_reslice_trans(:,:,z),-rotAngleLo,'bicubic','crop');
end
Y_F_lo_reslice_trans_rot(Y_F_lo_reslice_trans_rot<0) = 0.0;
Y_F_lo_reslice_trans_rot(Y_F_lo_reslice_trans_rot>niiFloMax) = niiFloMax;

% for z = 1:120
%     Y_F_lo_reslice_trans_rot(:,:,z) = imrotate(Y_F_lo_reslice_trans(:,:,z),-rotAngleLo,'bilinear','crop')';
% end
% 
% temp_dim = V_F_lo_reslice.dim;
% V_F_lo_reslice.dim(1) = temp_dim(2);
% V_F_lo_reslice.dim(2) = temp_dim(1);

%figure, imshow(imrotate(flipud(Y_F_lo_reslice_trans_rot(:,:,showSliceLo)),90),[]);

% first_lower_l - scanOffset
% figure, imshow(flipud(Y_F_up_reslice_trans_rot(:,:,first_lower_l - scanOffset)),[]);
% first_lower_l
% figure, imshow(flipud(Y_F_lo_reslice_trans_rot(:,:,first_lower_l)),[]);

V_F_lo_reslice.mat(1,4) = sign(V_F_lo_reslice.mat(1,4))*(abs(V_F_lo_reslice.mat(1,4)) + translationLo(1)*V_F_lo_reslice.mat(1,1));
V_F_lo_reslice.mat(2,4) = sign(V_F_lo_reslice.mat(2,4))*(abs(V_F_lo_reslice.mat(2,4)) + translationLo(1)*V_F_lo_reslice.mat(2,2));

V_F_lo_reslice = spm_write_vol(V_F_lo_reslice,Y_F_lo_reslice_trans_rot);
V_F_lo_reslice.dim

% plane_angle_Lo = 0.25*(angleL_Lo + angleR_Lo + angleXZ_Lo)
% 
% scanLo2 = scanLo;
% scanLo2.voxel.size = [voxelSize_lo(2) voxelSize_lo(1) voxelSize_lo(3)];
% scanLo2.voxel.dimensions = [size(niiFup.img,2) size(niiFup.img,1) size(niiFup.img,3)];   
% scanLo2.tiltAngle = plane_angle_Lo;
% 
% scanLo2.voxel.dimensions
% [V_F_lo_reslice,Y_F_lo_reslice] = STANCE_reslice_volume(V_F_lo_reslice,scanLo2,[],true,true);


%figure, imshow(flipud(Y_F_lo_reslice(:,:,showSliceLo)),[]);

% %first_lower_l - scanOffset
% figure, imshow(flipud(Y_F_up_reslice(:,:,first_lower_l - scanOffset)),[]);
% %first_lower_l
% figure, imshow(flipud(Y_F_lo_reslice(:,:,first_lower_l)),[]);
% 
% for z = 1:120
%     Y_F_lo_reslice_trans(:,:,z) = Y_F_lo_reslice(:,:,z)';
% end
% V_F_lo_reslice.dim = temp_dim;
% 
% V_F_lo_reslice = spm_write_vol(V_F_lo_reslice,Y_F_lo_reslice_trans);

%figure, imshow(flipud(Y_F_lo_reslice_trans(:,:,showSliceLo)),[]);
%figure, imshow(flipud(Y_F_up_reslice_trans(:,:,showSliceUp)),[]);

% first_lower_l - scanOffset
% figure, imshow(flipud(Y_F_up_reslice_trans(:,:,first_lower_l - scanOffset)),[]);
% first_lower_l
%figure, imshow(flipud(Y_F_lo_reslice_trans(:,:,first_lower_l)),[]);

Blo_reslice = Make_DIXON_threshold_image(niiFlo_reslice.img,niiWlo_reslice.img);
[tiltAngleLo_reslice,rotAngleLo_reslice] = Find_table_tilt_angle(Blo_reslice, voxelSize_lo,[5,50],[],false); %[20,60] %[55,115]  %[20,first_lower_l]
tiltAngleLo_reslice
rotAngleLo_reslice

Y_W_lo_reslice_trans = imtranslate(Y_W_lo_reslice,translationLo,'cubic','OutputView','same');
for z = 1:120
    Y_W_lo_reslice_trans_rot(:,:,z) = imrotate(Y_W_lo_reslice_trans(:,:,z),-rotAngleLo,'bicubic','crop');
end
Y_W_lo_reslice_trans_rot(Y_W_lo_reslice_trans_rot<0) = 0.0;
Y_W_lo_reslice_trans_rot(Y_W_lo_reslice_trans_rot>niiWloMax) = niiWloMax;


V_W_lo_reslice.dim
V_W_lo_reslice.mat(1,4) = V_F_lo_reslice.mat(1,4);
V_W_lo_reslice.mat(2,4) = V_F_lo_reslice.mat(2,4);

V_W_lo_reslice = spm_write_vol(V_W_lo_reslice,Y_W_lo_reslice_trans_rot);

% 
% Y_W_lo_reslice_trans = imtranslate(Y_W_lo_reslice,translationLo,'linear','OutputView','same');
% for z = 1:120
%     Y_W_lo_reslice_trans_rot(:,:,z) = imrotate(Y_W_lo_reslice_trans(:,:,z),-rotAngleLo,'bilinear','crop')';
% end
% V_W_lo_reslice.dim(1) = temp_dim(2);
% V_W_lo_reslice.dim(2) = temp_dim(1);
% V_W_lo_reslice.dim
% V_W_lo_reslice.mat(1,4) = V_F_lo_reslice.mat(1,4);
% V_W_lo_reslice.mat(2,4) = V_F_lo_reslice.mat(2,4);
% 
% V_W_lo_reslice = spm_write_vol(V_W_lo_reslice,Y_W_lo_reslice_trans_rot);
%
% [V_W_lo_reslice,Y_W_lo_reslice] = STANCE_reslice_volume(V_W_lo_reslice,scanLo2,[],true,true);
% 
% for z = 1:120
%     Y_W_lo_reslice_trans(:,:,z) = Y_W_lo_reslice(:,:,z)';
% end
% V_W_lo_reslice.dim = temp_dim;
% V_W_lo_reslice = spm_write_vol(V_W_lo_reslice,Y_W_lo_reslice_trans);

% Blo_reslice = Make_DIXON_threshold_image(niiFlo_reslice.img,niiWlo_reslice.img);
% 
% [tiltAngleLo_reslice,rotAngleLo_reslice] = Find_table_tilt_angle(Blo_reslice, voxelSize_lo,[5,50],[],false); %[30,60] %[20,first_lower_l]
% tiltAngleLo_reslice
% rotAngleLo_reslice
% [x_cm_Lo, angleL_Lo, angleR_Lo] = Find_plane_angle(Blo_reslice, voxelSize_lo, 2.5,[55,115],[],true);
% x_cm_Lo
% angleL_Lo
% angleR_Lo
% 
% [x_cm_Los, y_cm_Los, angleXZ_Lo, angleYZ_Lo] = Find_cm_angle(Bup_reslice, voxelSize_up, 2.5,[5,71],[],true);
% angleXZ_Lo
% %angleYZ_Lo
% 
% plane_angle_Lo_reslice = 0.25*(angleL_Lo + angleR_Lo + angleXZ_Lo)

V_W_lo_reslice.mat

% load rectified data
niiWlo_reslice = load_nii(NIFTI_file_name_W_lower_reslice);
figure, imshow(fliplr(rot90(niiWlo_reslice.img(:,:,showSliceLo))),[]);

figure, imshow(fliplr(rot90(niiFlo.img(:,:,showSliceLo))),[]);

niiFlo_reslice = load_nii(NIFTI_file_name_F_lower_reslice);
figure, imshow(fliplr(rot90(niiFlo_reslice.img(:,:,showSliceLo))),[]);

figure, imshow(fliplr(rot90(niiFup.img(:,:,showSliceUp))),[]);

figure, imshow(fliplr(rot90(niiFup_reslice.img(:,:,showSliceUp))),[]);

%scanOffset = round((V_W_up_reslice.mat(3,4) - V_W_lo_reslice.mat(3,4))/2.5)

% first_lower_l - scanOffset
% figure, imshow(flipud(Y_F_up_reslice_trans_rot(:,:,first_lower_l - scanOffset)),[]);
% first_lower_l
% figure, imshow(flipud(Y_F_lo_reslice_trans_rot(:,:,first_lower_l)),[]);

%% clean up memory and write output

varlist = {'niiFup', 'niiFlo','niiWup', 'niiWlo','niiFup_reslice', 'niiFlo_reslice','niiWup_reslice', 'niiWlo_reslice','V_F_up_reslice','V_F_lo_reslice','V_W_up_reslice','V_W_lo_reslice','V_F_up','V_F_lo','V_W_up','V_W_lo','Y_F_up_reslice','Y_F_lo_reslice','Y_W_up_reslice','Y_W_lo_reslice','Y_F_up','Y_F_lo','Y_W_up','Y_W_lo','Bup','Blo','Bup_reslice','Blo_reslice','Y_F_up_reslice_trans','Y_W_up_reslice_trans','Y_F_lo_reslice_trans','Y_W_lo_reslice_trans','Y_F_up_reslice_trans_rot','Y_W_up_reslice_trans_rot','Y_F_lo_reslice_trans_rot','Y_W_lo_reslice_trans_rot'};
clear(varlist{:})

delete(NIFTI_file_name_F_upper_NII, NIFTI_file_name_W_upper_NII, NIFTI_file_name_F_lower_NII, NIFTI_file_name_W_lower_NII);

save subject02_f_registered.mat;