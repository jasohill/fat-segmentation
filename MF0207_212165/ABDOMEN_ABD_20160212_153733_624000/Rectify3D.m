%------------------------------------------------------------------------
% Author: Dr. Jason E. Hill, post-doctoral researcher
% Institution: CNG at TTU in a working association with TTNI
% Date: 9 MAR 2015
% Updated: 9 MAR 2016
%------------------------------------------------------------------------

clear all, close all, clc;

addpath(genpath(pwd))
addpath('../../auxillary')
addpath('../../auxillary/NIfTI_20140122')
addpath('C:/spm/spm8')

% subject: MF0207
% session: 212165
% date:    12 FEB 2016

% image generation control
showProgress = false; %true;
batchShow    = ~showProgress;

%% scan info

NIFTI_file_name_F_upper = '20160212_153733t1vibedixontrap4bh320s004a1001.nii.gz';
NIFTI_file_name_W_upper = '20160212_153733t1vibedixontrap4bh320s005a1001.nii.gz';
NIFTI_file_name_F_lower = '20160212_153733t1vibedixontrap4bh320s009a1001.nii.gz';
NIFTI_file_name_W_lower = '20160212_153733t1vibedixontrap4bh320s010a1001.nii.gz';

NIFTI_file_name_F_upper_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_F_0004/20160212_153733t1vibedixontrap4bh320s004a1001.nii';
NIFTI_file_name_W_upper_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_W_0005/20160212_153733t1vibedixontrap4bh320s005a1001.nii';
NIFTI_file_name_F_lower_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_F_0009/20160212_153733t1vibedixontrap4bh320s009a1001.nii';
NIFTI_file_name_W_lower_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_W_0010/20160212_153733t1vibedixontrap4bh320s010a1001.nii';

NIFTI_file_name_F_upper_reslice = '20160212_153733t1vibedixontrap4bh320s004a1001_resliced.nii';
NIFTI_file_name_W_upper_reslice = '20160212_153733t1vibedixontrap4bh320s005a1001_resliced.nii';
NIFTI_file_name_F_lower_reslice = '20160212_153733t1vibedixontrap4bh320s009a1001_resliced.nii';
NIFTI_file_name_W_lower_reslice = '20160212_153733t1vibedixontrap4bh320s010a1001_resliced.nii';

voxelSize_up = [1.40625, 1.40625, 2.5];
interslice_spacing_fraction = 0.2;
voxelVolume_up = prod(voxelSize_up)*(1+interslice_spacing_fraction);
voxelSize_lo = [1.40625, 1.40625, 2.5];
voxelVolume_lo = prod(voxelSize_lo)*(1+interslice_spacing_fraction);
%X_shift = -round(2.2/voxelSize_lower(1));
%Y_shift = -round(9.1/voxelSize_lower(2));

% matching slices for stitching purposes
showSliceUp = 9;
showSliceLo = 115;

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

%addpath(genpath('C:\spm\spm8'))

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

% repeat for water saturated data
V_W_up = STANCE_load_volume(NIFTI_file_name_W_upper);
tic
[V_W_up_reslice,Y_W_up_reslice] = STANCE_reslice_volume(V_W_up,scanUp,[],true,true,3);
toc

niiWup_reslice = load_nii(NIFTI_file_name_W_upper_reslice);
Bup_reslice = Make_DIXON_threshold_image(niiFup_reslice.img,niiWup_reslice.img);

BupCoM = center_of_mass(Bup_reslice>0.5)
BupShowSliceCoM = center_of_mass(Bup_reslice(:,:,showSliceUp)>0.5)

translationUp = [BupShowSliceCoM(1)-261/2.0,...
                 BupShowSliceCoM(2)-321/2.0,0]

V_F_up_reslice.mat

tic
Y_F_up_reslice_trans = imtranslate(Y_F_up_reslice,translationUp,'cubic','OutputView','same');
toc

tic
Y_W_up_reslice_trans = imtranslate(Y_W_up_reslice,translationUp,'cubic','OutputView','same');
toc

figure, imshow(fliplr(rot90(flipud(Y_F_up_reslice_trans(:,:,showSliceUp)))),[]);

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

V_W_up_reslice.mat

% load rectified data

niiWup_reslice = load_nii(NIFTI_file_name_W_upper_reslice);
figure, imshow(fliplr(rot90(niiWup_reslice.img(:,:,showSliceUp))),[]);

niiFup_reslice = load_nii(NIFTI_file_name_F_upper_reslice);
figure, imshow(fliplr(rot90(niiFup_reslice.img(:,:,showSliceUp))),[]);

%% register landmarks (key anatomical features)
sessionOffset  = 0;

upperBottom = 9;
Z_L3_L4u    = 1;
Z_L2_L3u    = 6;
Z_L1_L2u    = 18;  % middle of L1-L2 intervertebral disk slice (upper scan)
% Z_L1b       = 29;  % bottom of L1 vertebra slice
Z_T4_T5t    = 104; % top of T4-T5 intervertebral disk slice
Z_T9_T10t   = 64;  % top of T10 vertebra disk slice (use for bottom of thoracic cavity) 
Z_T10_T11   = 53; 
upperTop = 113;

% NOTE: diaphraim is at slices = 
% - register from the bottom of T8
heartMax  = upperTop - 70+1;
heartApex = upperTop - 55+1;
diaphram  = upperTop - Z_T10_T11+1; % used as bottom of lung/thorax

upperSlices = upperTop:-1:upperBottom;

lowerBottom = 5;
FH       = 20;  % widest slice of the formoral head
Z_L5_S1b = 57;  % bottom of L5-S1 intervertebral disk slice
Z_L4_L5  = 72;   % middle of L4-L5 intervertebral disk slice
Z_L3_L4  = 86;   % middle of L3-L4 intervertebral disk slice
UM       = 92;   % umbilicus (belly button)
Z_L2_L3  = 99;  % middle of L2-L3 intervertebral disk slice
Z_L1_L2  = 111;  % middle of L1-L2 intervertebral disk slice
%Z_L1_L2t    = 118; % top of L2 vertebra disk
lowerTop = 117;
lastSlice = 36;  % lowest slice without clipping

lowerSlices = lowerTop:-1:lowerBottom;

% correct scan offset
scanOffset = 104   %Z_L1_L1 - Z_L1_L1u

first_lower_l = 121 - lowerTop + upperBottom + scanOffset;

Nslices = length(upperSlices) + length(lowerSlices);% + repeatEstimators;
slices = [upperSlices, lowerSlices];
levels(1:length(upperSlices))         = 1;
levels(length(upperSlices)+1:Nslices) = 2;
voxelVolumes(1:length(upperSlices))         = voxelVolume_up;
voxelVolumes(length(upperSlices)+1:Nslices) = voxelVolume_lo;

% NOTE: belly button is from i = 110-117
bellybutton = (length(upperSlices) + lowerTop - UM + 1 - 6):(length(upperSlices) + lowerTop - UM + 1 + 5);
bbx = 78;

% parameter settings
cutLineLength = 92;

backgroundThreshold   = 10; %15%10
vetoFactorF = 4.5;
vetoFactorW = 3.5;
foregroundThreshold1  = 45; %74; %95; %84; %74; %59
foregroundThreshold2  = 10;
foregroundThresholdF2 = 145;
foregroundThresholdV  = 128; %185; %175; %103;
foregroundThresholdV2 = 230;
foregroundThresholdW = 40;

nbrThresholdsF = 4; %4; %7;
nbrThresholdsW = 7; %4; %9;
nbrThresholdsV = 7; %4; %7;

seSCAT = strel('rectangle',[7,17]);  
VATprelim      = 17;     % opening for preliminary VAT segmentation

dims = size(niiFup.img);
% correction of inner mask by removing effect of mammaries
innerVeto = zeros(1,dims(3));
innerVeto(115:120) = 261 - 177;
innerVeto(114)     = 261 - 179;
innerVeto(113)     = 261 - 180;
innerVeto(112)     = 261 - 182;
innerVeto(111)     = 261 - 184;
innerVeto(110)     = 261 - 186;
innerVeto(109)     = 261 - 188;
innerVeto(108)     = 261 - 189;
innerVeto(107)     = 261 - 190;
innerVeto(106)     = 261 - 191;
innerVeto(105)     = 261 - 192;
innerVeto(104)     = 261 - 193;
innerVeto(103)     = 261 - 194;
innerVeto(102)     = 261 - 195;
innerVeto(101)     = 261 - 196;
innerVeto(100)     = 261 - 197;
innerVeto(99)      = 261 - 198;
innerVeto(98)      = 261 - 199;
innerVeto(97)      = 261 - 200;
innerVeto(96)      = 261 - 201;
innerVeto(95)      = 261 - 202;
innerVeto(94)      = 261 - 203;
innerVeto(93)      = 261 - 204;
innerVeto(92)      = 261 - 205;
innerVeto(88:91)   = 261 - 206;
innerVeto(86:87)   = 261 - 207;
innerVeto(84:85)   = 261 - 208;
innerVeto(82:83)   = 261 - 209;
innerVeto(77:81)   = 261 - 210;
innerVeto(60:76)   = 261 - 211;
%innerVeto      = 66;     % used for correcting the body cavity (inner) mask
%innerVetoLimit = 22;
voidThreshold  = 5000;   % threshold of count of number of void pixels
lungErode      = 6; %6

aortaThreshold = 80;     % threshold of 8-bit water dominated signal data
aortaMinArea   = 70;
aortaDistance  = 50;
aortaRadius    = 9;

notLungDilate  = 35; %45%42%21;
innerDilate    = 10; %13;
heartShift     = 30;
heartDilate    = 45; %45 %42
heartAreaThreshold = 400;
CATareaThreshold = 64; %160;

CATdilate      = 13; %7; %5; %3; %4; %6
CATmargin      = 3; %6;
PAATdilate     = 10; %6;

% CATdilateCorrect = zeros(1,Nslices);
% CATdilateCorrect(heartMax-24:heartMax-19) = +3; 

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
aortaSeed(261-126,164,40:43) = 1;
aortaSeed(261-124,164,44:47) = 1;
aortaSeed(261-123,163,48:51) = 1;
aortaSeed(261-122,163,52:54) = 1;
aortaSeed(261-121,163,55:58) = 1;
aortaSeed(261-120,163,58:59) = 1;
aortaSeed(261-119,164,60:66) = 1;
aortaSeed(261-118,164,67:69) = 1;
aortaSeed(261-117,164,69)    = 1;
aortaSeed(261-117,164,70:71) = 1;
aortaSeed(261-115,164,72:73) = 1;
aortaSeed(261-116,164,74:76) = 1;
aortaSeed(261-115,165,77:80) = 1;
aortaSeed(261-114,164,81:84) = 1;
aortaSeed(261-114,165,85:100) = 1; 
aortaSeed(261-114,166,101:104)     = 1;
% curved section
aortaSeed(261-114,165,105:upperTop) = 1;
aortaSeed(261-116,165,104:upperTop) = 1;
aortaSeed(261-119,165,107:upperTop) = 1;
aortaSeed(261-130,165,107:upperTop) = 1;
aortaSeed(261-115,161,107:upperTop) = 1;
aortaSeed(261-139,156,107:upperTop) = 1;
aortaSeed(261-139,151,107:upperTop) = 1;
aortaSeed(261-153,149,107:upperTop) = 1;
aortaSeed(261-130,162,108:upperTop) = 1;

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

[tiltAngleLo,rotAngleLo] = Find_table_tilt_angle(Blo, voxelSize_lo,[5,65],[],false); %[20,60] %60 % [20,first_lower_l]
tiltAngleLo
rotAngleLo

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
Blo_reslice = Make_DIXON_threshold_image(niiFlo_reslice.img,niiWlo_reslice.img);
%figure, imshow(fliplr(Blo_reslice(:,:,showSliceLo)),[]);

BloCoM = center_of_mass(Blo_reslice>0.5)
BloShowSliceCoM = center_of_mass(Blo_reslice(:,:,showSliceLo)>0.5)

translationLo = [BloShowSliceCoM(1)-261/2.0,...
                 BloShowSliceCoM(2)-321/2.0,0]
% NOTE: translation = [50, 0, 0]; shifts up
% NOTE: translation = [0, 50, 0]; shifts right

Y_F_lo_reslice_trans = imtranslate(Y_F_lo_reslice,translationLo,'cubic','OutputView','same');

for z = 1:120
    Y_F_lo_reslice_trans_rot(:,:,z) = imrotate(Y_F_lo_reslice_trans(:,:,z),-rotAngleLo,'bicubic','crop');
end
Y_F_lo_reslice_trans_rot(Y_F_lo_reslice_trans_rot<0) = 0.0;
Y_F_lo_reslice_trans_rot(Y_F_lo_reslice_trans_rot>niiFloMax) = niiFloMax;

V_F_lo_reslice.mat(1,4) = sign(V_F_lo_reslice.mat(1,4))*(abs(V_F_lo_reslice.mat(1,4)) + translationLo(1)*V_F_lo_reslice.mat(1,1));
V_F_lo_reslice.mat(2,4) = sign(V_F_lo_reslice.mat(2,4))*(abs(V_F_lo_reslice.mat(2,4)) + translationLo(2)*V_F_lo_reslice.mat(2,2));

V_F_lo_reslice = spm_write_vol(V_F_lo_reslice,Y_F_lo_reslice_trans_rot);

Blo_reslice = Make_DIXON_threshold_image(niiFlo_reslice.img,niiWlo_reslice.img);
[tiltAngleLo_reslice,rotAngleLo_reslice] = Find_table_tilt_angle(Blo_reslice, voxelSize_lo,[5,65],[],false); %[20,60] %[55,115]  %[20,first_lower_l]
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

V_W_lo_reslice.mat

% load rectified data
niiWlo_reslice = load_nii(NIFTI_file_name_W_lower_reslice);
figure, imshow(fliplr(rot90(niiWlo_reslice.img(:,:,showSliceLo))),[]);

figure, imshow(fliplr(rot90(niiFlo.img(:,:,showSliceLo))),[]);

niiFlo_reslice = load_nii(NIFTI_file_name_F_lower_reslice);
figure, imshow(fliplr(rot90(niiFlo_reslice.img(:,:,showSliceLo))),[]);

figure, imshow(fliplr(rot90(niiFup.img(:,:,showSliceUp))),[]);

figure, imshow(fliplr(rot90(niiFup_reslice.img(:,:,showSliceUp))),[]);

%% clean up memory and write output

varlist = {'niiFup', 'niiFlo','niiWup', 'niiWlo','niiFup_reslice', 'niiFlo_reslice','niiWup_reslice', 'niiWlo_reslice','V_F_up_reslice','V_F_lo_reslice','V_W_up_reslice','V_W_lo_reslice','V_F_up','V_F_lo','V_W_up','V_W_lo','Y_F_up_reslice','Y_F_lo_reslice','Y_W_up_reslice','Y_W_lo_reslice','Y_F_up','Y_F_lo','Y_W_up','Y_W_lo','Bup','Blo','Bup_reslice','Blo_reslice','Y_F_up_reslice_trans','Y_W_up_reslice_trans','Y_F_lo_reslice_trans','Y_W_lo_reslice_trans','Y_F_up_reslice_trans_rot','Y_W_up_reslice_trans_rot','Y_F_lo_reslice_trans_rot','Y_W_lo_reslice_trans_rot'};
clear(varlist{:})

delete(NIFTI_file_name_F_upper_NII, NIFTI_file_name_W_upper_NII, NIFTI_file_name_F_lower_NII, NIFTI_file_name_W_lower_NII);

save subject07_i_registered.mat;