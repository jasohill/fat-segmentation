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
% session: 212165_2
% date:    19 MAR 2016

%% scan info

NIFTI_file_name_F_upper = '20160310_135937t1vibedixontrap4bh320s004a1001.nii.gz';
NIFTI_file_name_W_upper = '20160310_135937t1vibedixontrap4bh320s005a1001.nii.gz';
NIFTI_file_name_F_lower = '20160310_135937t1vibedixontrap4bh320s009a1001.nii.gz';
NIFTI_file_name_W_lower = '20160310_135937t1vibedixontrap4bh320s010a1001.nii.gz';

NIFTI_file_name_F_upper_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_F_0004/20160310_135937t1vibedixontrap4bh320s004a1001.nii';
NIFTI_file_name_W_upper_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_W_0005/20160310_135937t1vibedixontrap4bh320s005a1001.nii';
NIFTI_file_name_F_lower_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_F_0009/20160310_135937t1vibedixontrap4bh320s009a1001.nii';
NIFTI_file_name_W_lower_NII = 'T1_VIBE_DIXON_TRA_P4_BH_320_W_0010/20160310_135937t1vibedixontrap4bh320s010a1001.nii';

NIFTI_file_name_F_upper_reslice = '20160310_135937t1vibedixontrap4bh320s004a1001_resliced.nii';
NIFTI_file_name_W_upper_reslice = '20160310_135937t1vibedixontrap4bh320s005a1001_resliced.nii';
NIFTI_file_name_F_lower_reslice = '20160310_135937t1vibedixontrap4bh320s009a1001_resliced.nii';
NIFTI_file_name_W_lower_reslice = '20160310_135937t1vibedixontrap4bh320s010a1001_resliced.nii';

voxelSize_up = [1.40625, 1.40625, 2.5];
interslice_spacing_fraction = 0.2;
voxelVolume_up = prod(voxelSize_up)*(1+interslice_spacing_fraction);
voxelSize_lo = [1.40625, 1.40625, 2.5];
voxelVolume_lo = prod(voxelSize_lo)*(1+interslice_spacing_fraction);

% matching slices for stitching purposes
showSliceUp = 15;
showSliceLo = 105;

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
sessionOffset  = 4;

upperBottom = 10;
Z_L2_L3u    = 10; %L3 is quarter vertabra
Z_L1_L2u    = 22;  % middle of L1-L2 intervertebral disk slice (upper scan)
% Z_L1b        = 33;  % bottom of L1 vertebra slice
Z_T4_T5t    = 108; % top of T4-T5 intervertebral disk slice
%upperSlices = Z_T4_T5t:-1:Z_L1_L2u;
Z_T9_T10t      = 67;  % top of T10 vertebra disk slice (use for bottom of thoracic cavity) 
Z_T10_T11   = 57; 
upperTop = 110;

% NOTE: diaphraim is at slices = 
% - register from the bottom of T8
heartMax  = upperTop-73+1;
heartApex = upperTop-58+1;
diaphram  = upperTop-Z_T10_T11+1; % used as bottom of lung/thorax

upperSlices = upperTop:-1:upperBottom;

lowerBottom = 6;
FH       = 22;  % widest slice of the formoral head
Z_L5_S1b = 59;  % bottom of L5-S1 intervertebral disk slice
Z_L4_L5  = 73;   % middle of L4-L5 intervertebral disk slice
Z_L3_L4  = 87;   % middle of L3-L4 intervertebral disk slice
UM       = 92;   % umbilicus (belly button)
Z_L2_L3  = 100;  % middle of L2-L3 intervertebral disk slice %L3 is quarter vertabra
Z_L1_L2  = 112;  % middle of L1-L2 intervertebral disk slice
%Z_L1_L2t    = 118; % top of L2 vertebra disk
% Z_L1_T12 = 34; %X_L1_T12 = 163; Y_L1_T12 = 128;
% Z_T12_T11  = 46; %X_T12_T11 = 162; Y_T12_T11 = 104;
% Z_T11_T10  = 56; %X_T11_T10 = 162; Y_T11_T10 = 100;
% Z_T10_T9   = 63; %X_T10_T9 = 163; Y_T10_T9 = 98; %Mid Pt
% Z_T9_T8  = 76; %X_T9_T8 = 163; Y_T9_T8 = 94;
% Z_T8_T7  = 84; %X_T8_T7 = 163; Y_T8_T7 = 94;
% Z_T7_T6  = 93; %X_T7_T6 = 163; Y_T7_T6 = 90;
% Z_T6_T5  = 101; %X_T6_T5 = 164; Y_T6_T5 = 89; %T5 is quarter vertabra
% Z_T5_T4  = 108; %X_T5_T4 = 164; Y_T5_T4 = 91;
% Z_T4_T3  = 117; %X_T4_T3 = 164; Y_T4_T3 = 96;
% Z_T3_T2  =
% Z_T2_T1  = 
% Z_T1_C7  =
% Z_C7_C6  =
% Z_C6_C5  =
% Z_C5_C4  =
% Z_C4_C3  =
% Z_C3_C2  =
% Z_C2_C1  =
lowerTop = 116;
lastSlice = 57;  % lowest slice without clipping

lowerSlices = lowerTop:-1:lowerBottom;

% correct scan offset
scanOffset = 90   %Z_L1_L1 - Z_L1_L1u

first_lower_l = 121 - lowerTop + upperBottom + scanOffset;

Nslices = length(upperSlices) + length(lowerSlices);% + repeatEstimators;
slices = [upperSlices, lowerSlices];
levels(1:length(upperSlices))         = 1;
levels(length(upperSlices)+1:Nslices) = 2;
voxelVolumes(1:length(upperSlices))         = voxelVolume_up;
voxelVolumes(length(upperSlices)+1:Nslices) = voxelVolume_lo;

% NOTE: belly button is from i = 110-117
bellybutton = (length(upperSlices) + lowerTop - UM + 1 - 6):(length(upperSlices) + lowerTop - UM + 1 + 5);
bbx = 76;

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
foregroundThresholdW  = 40;

nbrThresholdsF = 4; %4; %7;
nbrThresholdsW = 7; %4; %9;
nbrThresholdsV = 7; %4; %7;

seSCAT = strel('rectangle',[7,17]);  
VATprelim      = 17;     % opening for preliminary VAT segmentation

dims = size(niiFup.img);
% correction of inner mask by removing effect of mammaries
innerVeto = zeros(1,dims(3));
innerVeto(114:120) = 261 - 186;
innerVeto(113)     = 261 - 188;
innerVeto(112)     = 261 - 189;
innerVeto(111)     = 261 - 190;
innerVeto(110)     = 261 - 191;
innerVeto(109)     = 261 - 192;
innerVeto(108)     = 261 - 193;
innerVeto(107)     = 261 - 194;
innerVeto(106)     = 261 - 195;
innerVeto(104:105) = 261 - 196;
innerVeto(103)     = 261 - 197;
innerVeto(102)     = 261 - 199;
innerVeto(101)     = 261 - 200;
innerVeto(100)     = 261 - 201;
innerVeto(99)      = 261 - 202;
innerVeto(98)      = 261 - 203;
innerVeto(97)      = 261 - 204;
innerVeto(96)      = 261 - 205;
innerVeto(95)      = 261 - 206;
innerVeto(92:94)   = 261 - 207;
innerVeto(88:91)   = 261 - 208;
innerVeto(86:87)   = 261 - 209;
innerVeto(84:85)   = 261 - 210;
innerVeto(82:83)   = 261 - 211;
innerVeto(77:81)   = 261 - 212;
innerVeto(60:76)   = 261 - 213;
%innerVeto      = 60;     % used for correcting the body cavity (inner) mask
%innerVetoLimit = 20;
voidThreshold  = 5000;   % threshold of count of number of void pixels
lungErode      = 10; %6; %6

aortaThreshold = 80;     % threshold of 8-bit water dominated signal data
aortaMinArea   = 70;
aortaDistance  = 50;
aortaRadius    = 9;

notLungDilate  = 40; %45%42%21;
innerDilate    = 11; %13;
heartShift     = 25;
heartDilate    = 45; %45 %42
heartAreaThreshold = 400;
CATareaThreshold = 64% 160;

CATdilate      = 13; %7; %5; %3; %4; %6
CATmargin      = 3; %6;
PAATdilate     = 6;

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

dims = size(niiFup.img);

% define aorta seed in native coordinates of level 1 (upper scan)
aortaSeed = false(dims(2),dims(1),dims(3));
% vertical section
aortaSeed(261-130,167,40:43)   = 1;
aortaSeed(261-129,165,44:46)   = 1;
aortaSeed(261-128,165,47:48)   = 1;
aortaSeed(261-127,166,49:51)   = 1;
aortaSeed(261-126,167,52:54)   = 1;
aortaSeed(261-125,167,55:58)   = 1;
aortaSeed(261-124,167,58:59)   = 1;
aortaSeed(261-121,169,60:66)   = 1;
aortaSeed(261-121,168,67:69)   = 1;
aortaSeed(261-121,169,70:71)   = 1;
aortaSeed(261-120,169,72:73)   = 1;
aortaSeed(261-119,169,74:76)   = 1;
aortaSeed(261-118,170,77:80)   = 1;
aortaSeed(261-117,171,81:82)   = 1;
aortaSeed(261-116,171,83:87)   = 1;
aortaSeed(261-115,172,88:100)   = 1; 
aortaSeed(261-116,172,101:104)     = 1;
% curved section
aortaSeed(261-118,171,105:upperTop) = 1;
aortaSeed(261-120,171,106:upperTop) = 1;
aortaSeed(261-123,171,107:upperTop) = 1;
aortaSeed(261-145,169,107:upperTop) = 1;
aortaSeed(261-150,167,107:upperTop) = 1;
aortaSeed(261-143,162,107:upperTop) = 1;
aortaSeed(261-143,157,107:upperTop) = 1;
aortaSeed(261-157,155,107:upperTop) = 1;
aortaSeed(261-134,168,108:upperTop) = 1;

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

[tiltAngleLo,rotAngleLo] = Find_table_tilt_angle(Blo, voxelSize_lo,[5,55],[],false); %[20,60] %60 % [20,first_lower_l]
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
[tiltAngleLo_reslice,rotAngleLo_reslice] = Find_table_tilt_angle(Blo_reslice, voxelSize_lo,[5,55],[],false); %[20,60] %[55,115]  %[20,first_lower_l]
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

save subject07_f_registered.mat;