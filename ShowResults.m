clear all; close all;
addpath(genpath(pwd));

subject = 15;
pre = true;
%pre = false;
pre_AND_post = true;
%pre_AND_post = false;
write_results = false; %true;
%render = 'none';
render = 'PAT';
showSlice = true;

if pre
   letter = 'i'; %#ok<UNRCH>
   letterOther = 'f';
   sheet = 1;
else
   letter = 'f';
   letterOther = 'i';
   sheet = 2;
end
if ~pre_AND_post
   letterOther = letter;
end
if subject < 10
    results_fn = ['subject0',num2str(subject),'_',letter,'_results.mat'];
    results_fn_other = ['subject0',num2str(subject),'_',letterOther,'_results.mat'];
else
    results_fn = ['subject',num2str(subject),'_',letter,'_results.mat'];
    results_fn_other = ['subject',num2str(subject),'_',letterOther,'_results.mat'];
end

% other scan for comparison
load(results_fn_other)
sessionOffsetOther = sessionOffset;
upperTopOther =  upperTop;
upperBottomOther =  upperBottom;
upperSlicesOther = upperSlices;
lowerTopOther =  lowerTop;
lowerBottomOther =  lowerBottom;
scanOffsetOther = scanOffset;
diaphramOther = diaphram;
lastSliceOther = lastSlice;
FHother = FH;

% scan to display:
load(results_fn)
sliceNbr = 46; %heartMax;

% spreadsheet options
rows = [3 4 5 NaN 6 NaN 7 8 9 10 11 12 13 14 15];
row = rows(subject);
if write_results
filename = 'volumes.xlsx';

if sheet == 1
A = {'Quantity','Number of overlapping slices','Number of Thorax slices','Number of abdomenal slices','Number of abdomenal slices no clip','Number of clipped slices','SCAT volume in thorax (mL)','SCAT volume in abdomen no clipping (mL)','SCAT volume in abdomen to femoral head (mL)','SCAT estimated volume in abdomen to femoral head (mL)','VAT plus volume in thorax (mL)','VAT plus volume in abdomen no clipping (mL)','VAT plus volume in abdomen to femoral head (mL)','Organs plus volume in thorax (mL)','Organs plus volume in abdomen no clipping (mL)','Organs plus volume in abdomen to femoral head (mL)','Voids volume in thorax (mL)','Voids volume in abdomen no clipping (mL)','Voids volume in abdomen to femoral head (mL)','Lung volume (mL)','Heart volume (mL)','CAT volume (mL)','Aorta volume (mL)','PAAT volume (mL)','EAT volume (mL)','PAT volume (mL)','Organ Average Fat Fraction','TAT volume (mL)','true VAT volume in abdomen no clip (mL)','true VAT volume to femoral head (mL)','IMAT volume in thorax (mL)','IMAT volume in abdomen no clip (mL)','IMAT volume in abdomen to femoral head (mL)';...
     'ID','NOpre','Nthorax','Nabdomen','NabdomenNC','Nclipped','SCATVTpre','SCATVANCpre','SCATVAFHpre','SCATestVAFHpre','VATPVTpre','VATPVANCpre','VATPVAFHpre','OPVTpre','OPVANCpre','OPVAFHpre','VVTpre','VVANCpre','VVAFHpre','LUNGVpre','HEARTVpre','CATVpre','AORTAVpre','PAATVpre','EATVpre','PATVpre','OaveFFpre','TATVpre','VATVANCpre','VATVFHpre','IMATVTpre','IMATVANCpre','IMATVAFHpre'};
elseif sheet == 2
A = {'Quantity','Number of overlapping slices','Number of Thorax slices','Number of abdomenal slices','Number of abdomenal slices no clip','Number of clipped slices','SCAT volume in thorax (mL)','SCAT volume in abdomen no clipping (mL)','SCAT volume in abdomen to femoral head (mL)','SCAT estimated volume in abdomen to femoral head (mL)','VAT plus volume in thorax (mL)','VAT plus volume in abdomen no clipping (mL)','VAT plus volume in abdomen to femoral head (mL)','Organs plus volume in thorax (mL)','Organs plus volume in abdomen no clipping (mL)','Organs plus volume in abdomen to femoral head (mL)','Voids volume in thorax (mL)','Voids volume in abdomen no clipping (mL)','Voids volume in abdomen to femoral head (mL)','Lung volume (mL)','Heart volume (mL)','CAT volume (mL)','Aorta volume (mL)','PAAT volume (mL)','EAT volume (mL)','PAT volume (mL)','Organ Average Fat Fraction','TAT volume (mL)','true VAT volume in abdomen no clip (mL)','true VAT volume to femoral head (mL)','IMAT volume in thorax (mL)','IMAT volume in abdomen no clip (mL)','IMAT volume in abdomen to femoral head (mL)';...
     'ID','NOpost','Nthorax','Nabdomen','NabdomenNC','Nclipped','SCATVTpost','SCATVANCpost','SCATVAFHpost','SCATestVAFHpost','VATPVTpost','VATPVANCpost','VATPVAFHpost','OPVTpost','OPVANCpost','OPVAFHpost','VVTpost','VVANCpost','VVAFHpost','LUNGVpost','HEARTVpost','CATVpost','AORTAVpost','PAATVpost','EATVpost','PATVpost','OaveFFpost','TATVpost','VATVANCpost','VATVFHpost','IMATVTpost','IMATVANCpost','IMATVAFHpost'};
end
xlRange = 'A1';
xlswrite(filename,A,sheet,xlRange)

A = {1;2;3;5;7;8;9;10;11;12;13;14;15};
xlRange = 'A3';
xlswrite(filename,A,sheet,xlRange)
end

%render = 'none';
%render = 'VAT';
%render = 'organs';
%render = 'lung';
%render = 'heart';
%render = 'CAT';
%render = 'aorta';
%render = 'PAAT';

%lowerTopLevel1 = 121 + scanOffset - lowerTop;
%lowerTopLevel1 = 90 + upperBottom - 1;
%lastSlice  = 58;

%first_lower_i = length(upperSlices) + 122 - upperBottom - lowerTopLevel1;
%first_lower_i = 2*length(upperSlices) - (lowerTopLevel1 - upperBottom + 1) + 1;
%first_lower_i = 2*length(upperSlices) - (lowerTopLevel1) + 1;
first_lower_i = length(upperSlices) + lowerTop - upperBottom - scanOffset;

show_i = [1:length(upperSlices) first_lower_i:Nslices];
%diaphram_i = [diaphram:length(upperSlices) first_lower_i:Nslices];
overlap_lower_i = (length(upperSlices)+1):(first_lower_i-1);
overlap_upper_i = (length(upperSlices)-length(overlap_lower_i)+1):length(upperSlices);
no_overlap_i = [1:(length(upperSlices)-length(overlap_lower_i)) first_lower_i:Nslices];

Nshow = length(show_i)

Noverlap = length(overlap_lower_i)
% reconstruct TATvolume
for i = 1:length(TATvolume)
    voxelVolume = voxelVolumes(i);    
    TATmask = (TAT3d(:,:,i)>0);
    TATvolume(i)   = double(sum(sum(TATmask)))*voxelVolume/1000.0;   % [mL]
end
TATvolume(length(TATvolume)+1:length(SCATvolume)) = 0;
% take a weighted average of volumes from both scans for overlaping region
for j = 1:Noverlap
    average_weight = j/(Noverlap+1);
    SCATvolume(overlap_upper_i(j))   = average_weight*SCATvolume(overlap_lower_i(j))   + (1-average_weight)*SCATvolume(overlap_upper_i(j));
    VATvolume(overlap_upper_i(j))    = average_weight*VATvolume(overlap_lower_i(j))    + (1-average_weight)*VATvolume(overlap_upper_i(j));
    IMATvolume(overlap_upper_i(j))   = average_weight*IMATvolume(overlap_lower_i(j))   + (1-average_weight)*IMATvolume(overlap_upper_i(j));
    organsVolume(overlap_upper_i(j)) = average_weight*organsVolume(overlap_lower_i(j)) + (1-average_weight)*organsVolume(overlap_upper_i(j));
    voidsVolume(overlap_upper_i(j))  = average_weight*voidsVolume(overlap_lower_i(j))  + (1-average_weight)*voidsVolume(overlap_upper_i(j));       
end

% compute shared thorax slices
Nthorax = min(diaphram,diaphramOther)
if diaphram > diaphramOther
    thorax_i = [(diaphram-Nthorax+1):diaphram];

else
    thorax_i = [1:diaphram];
    sliceNbr = sliceNbr - (diaphramOther-Nthorax);
end

% compute shared abdomen slices (with clipping allowed)
FH_i = length(upperSlices)+lowerTop-FH+1
diaphram 
NabdomenThis = FH_i - diaphram - Noverlap

FHother_i = length(upperSlicesOther)+lowerTopOther-FHother+1
NoverlapOther = max(lowerTopOther - upperBottomOther - scanOffsetOther - 1,0)
NabdomenOther = FHother_i - diaphramOther - NoverlapOther

Nabdomen = min(NabdomenThis,NabdomenOther)
if NabdomenThis > NabdomenOther
    abdomen_i = [(diaphram+2):length(upperSlices) first_lower_i:(FH_i-(NabdomenThis-NabdomenOther-1))];
else
    abdomen_i = [(diaphram+1):length(upperSlices) first_lower_i:FH_i];
end

% handle clipped slices whose areas fall outside of the FOV
if lastSlice>0
    Nclipped = max(lastSlice - FH,lastSliceOther - FHother)
    abdomen_no_clip_i =  abdomen_i(1:end-Nclipped);
    abdomen_clipped_i =  abdomen_i(end-Nclipped+1:end);
end
NabdomenNC = length(abdomen_no_clip_i)

if write_results
A = [Noverlap Nthorax Nabdomen NabdomenNC Nclipped];
xlRange = strcat('B',num2str(row));
xlswrite(filename,A,sheet,xlRange)
end

% show results
figure, 
imshow(Isaggital(show_i,:),[])
daspect([3 1.40625 1])
%end_i = length(SCATvolume)-repeatEstimators;

EATvolume(Nthorax+1:length(CATvolume)) = 0;
yT = [aortaVolume(show_i);  PAATvolume(show_i); heartVolume(show_i); EATvolume(show_i); CATvolume(show_i); voidsVolume(show_i); lungVolume(show_i); TATvolume(show_i); VATvolume(show_i); organsVolume(show_i); IMATvolume(show_i); SCATvolume(show_i)]'; 
figure, 
bhT = barh(yT,'stacked','LineStyle','none');
set(bhT(1),'FaceColor','magenta');
set(bhT(2),'FaceColor','cyan');
set(bhT(3),'FaceColor','red');
set(bhT(4),'FaceColor',[0.75,0.75,0.75]);
set(bhT(5),'FaceColor',[0,0.75,0.5]);
set(bhT(6),'FaceColor',[0,0,0.5]);
set(bhT(7),'FaceColor','blue');
set(bhT(8),'FaceColor',[0,0.5,0]);
set(bhT(9),'FaceColor','green');
set(bhT(10),'FaceColor',[0.5,0,0]);
set(bhT(11),'FaceColor',[1 0.667 0.333]);
set(bhT(12),'FaceColor','yellow');
set(gca,'Color','k')
axT = gca; 
axT.YDir = 'reverse';
axT.YLim = [1 Nshow];
xlabel('Volume [mL]')
ylabel('slice number')
title('Total tissue segmentation')

yFp = [PAATvolume(show_i); EATvolume(show_i); CATvolume(show_i); TATvolume(show_i) + VATvolume(show_i) + IMATvolume(show_i); SCATvolume(show_i)]'; 
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

yF = [PAATvolume(show_i); EATvolume(show_i); CATvolume(show_i); TATvolume(show_i); VATvolume(show_i); IMATvolume(show_i); SCATvolume(show_i)]'; 
figure, 
bhF = barh(yF,'stacked','LineStyle','none');
set(bhF(1),'FaceColor','cyan');
set(bhF(2),'FaceColor',[0.75,0.75,0.75]);
set(bhF(3),'FaceColor',[0,0.75,0.5]);
set(bhF(4),'FaceColor',[0,0.5,0]);
set(bhF(5),'FaceColor','green');
set(bhF(6),'FaceColor',[1 0.667 0.333]);
set(bhF(7),'FaceColor','yellow');
set(gca,'Color','k')
axF = gca; 
axF.YDir = 'reverse';
axF.YLim = [1 Nshow];
xlabel('Volume [mL]')
ylabel('slice number')
title('Adipose tissue segmentation')

yFreport = [CATvolume(show_i); VATvolume(show_i); SCATvolume(show_i)]'; 
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


yW = [aortaVolume(show_i); heartVolume(show_i); organsVolume(show_i)]'; 
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

yV = [lungVolume(show_i); voidsVolume(show_i)]'; 
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

yCF = [PAATvolume(thorax_i); EATvolume(thorax_i); CATvolume(thorax_i)]'; 
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

yCV = [aortaVolume(thorax_i); PAATvolume(thorax_i); heartVolume(thorax_i); EATvolume(thorax_i); CATvolume(thorax_i)]'; 
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

disp('--------------')
disp('SCAT volumes: ')
disp('--------------')
SCATthorax = sum(SCATvolume(thorax_i));
X = ['Thoracic           : ',num2str(SCATthorax)];
disp(X)
SCATabNoClip = sum(SCATvolume(abdomen_no_clip_i));
X = ['Abdominal (no clip): ',num2str(SCATabNoClip)];
disp(X)
SCATabClipped = sum(SCATvolume(abdomen_i));
X = ['Abdominal to FH:     ',num2str(SCATabClipped)];
disp(X)
% estimate SCAT for clipped slices
for k = 1:Nclipped
    SCATbinary  = (SCAT3d(:,:,abdomen_clipped_i(k))>0);
    SCATsum1    = sum(SCATbinary);
    SCATsum2    = sum(fliplr(SCATbinary));
    SCATcorrect = sum((SCATsum1==0).*SCATsum2);
    SCATvolume(abdomen_clipped_i(k)) = SCATvolume(abdomen_clipped_i(k)) + SCATcorrect.*voxelVolumes(abdomen_clipped_i(k))/1000.0;
end
SCATest = sum(SCATvolume(abdomen_i));
X = ['Estimated abdominal to FH:     ',num2str(SCATest)];
disp(X)
if write_results
A = [SCATthorax SCATabNoClip SCATabClipped SCATest];
xlRange = strcat('G',num2str(row));
xlswrite(filename,A,sheet,xlRange)
end

disp('--------------')
disp('VAT+ volumes: ')
disp('--------------')
VATPthorax = sum(VATvolume(thorax_i)+IMATvolume(thorax_i));
X = ['Thoracic           : ',num2str(VATPthorax)];
disp(X)
VATPabNoClip = sum(VATvolume(abdomen_no_clip_i)+IMATvolume(abdomen_no_clip_i));
X = ['Abdominal(no clip): ',num2str(VATPabNoClip)];
disp(X)
VATPabdomen = sum(VATvolume(abdomen_i)+IMATvolume(abdomen_i));
X = ['Abdominalto FH:     ',num2str(VATPabdomen)];
disp(X)
if write_results
A = [VATPthorax VATPabNoClip VATPabdomen];
xlRange = strcat('K',num2str(row));
xlswrite(filename,A,sheet,xlRange)
end

disp('--------------')
disp('Organ+ volumes: ')
disp('--------------')
organsThorax = sum(organsVolume(thorax_i));
X = ['Thoracic           : ',num2str(organsThorax)];
disp(X)
organsAbNoClip = sum(organsVolume(abdomen_no_clip_i));
X = ['Abdominal(no clip): ',num2str(organsAbNoClip)];
disp(X)
organsAbdomen = sum(organsVolume(abdomen_i));
X = ['Abdominalto FH:     ',num2str(organsAbdomen )];
disp(X)
if write_results
A = [organsThorax organsAbNoClip organsAbdomen];
xlRange = strcat('N',num2str(row));
xlswrite(filename,A,sheet,xlRange)
end

disp('--------------')
disp('Voids volumes: ')
disp('--------------')
voidsThorax = sum(voidsVolume(thorax_i));
X = ['Thoracic           : ',num2str(voidsThorax)];
disp(X)
voidsAbNoClip = sum(voidsVolume(abdomen_no_clip_i));
X = ['Abdominal(no clip): ',num2str(voidsAbNoClip)];
disp(X)
voidsAbdomen = sum(voidsVolume(abdomen_i));
X = ['Abdominalto FH:     ',num2str(voidsAbdomen)];
disp(X)
if write_results
A = [voidsThorax voidsAbNoClip voidsAbdomen];
xlRange = strcat('Q',num2str(row));
xlswrite(filename,A,sheet,xlRange)
end

disp('--------------')
disp('Lung volume: ')
disp('--------------')
lungsVol = sum(lungVolume(thorax_i));
X = ['Thoracic           : ',num2str(lungsVol)];
disp(X)
disp('--------------')
disp('Heart volume: ')
disp('--------------')
heartVol = sum(heartVolume(thorax_i));
X = ['Thoracic           : ',num2str(heartVol)];
disp(X)
disp('--------------')
disp('CAT volume: ')
disp('--------------')
EATvol = sum(EATvolume(thorax_i));
PATvol = sum(CATvolume(thorax_i));
CATvol = EATvol + PATvol;
X = ['Thoracic           : ',num2str(CATvol)];
disp(X)
disp('--------------')
disp('Aorta volume: ')
disp('--------------')
aortaVol = sum(aortaVolume(thorax_i));
X = ['Thoracic           : ',num2str(aortaVol)];
disp(X)
disp('--------------')
disp('PAAT volume: ')
disp('--------------')
PAATvol = sum(PAATvolume(thorax_i));
X = ['Thoracic           : ',num2str(PAATvol)];
disp(X)
disp('--------------')
disp('EAT volume: ')
disp('--------------')
X = ['Thoracic           : ',num2str(EATvol)];
disp(X)
X = ['est error       +/-: ',num2str(2.0*abs(EATvol - sum(medfilt1(EATvolume(thorax_i)))))];
disp(X)
disp('--------------')
disp('PAT volume: ')
disp('--------------')
X = ['Thoracic           : ',num2str(PATvol)];
disp(X)
X = ['est error       +/-: ',num2str(2.0*abs(PATvol - sum(medfilt1(CATvolume(thorax_i)))))];
disp(X)
if write_results
A = [lungsVol, heartVol, CATvol, aortaVol, PAATvol, EATvol, PATvol];
xlRange = strcat('T',num2str(row));
xlswrite(filename,A,sheet,xlRange)
end

%FF_thorax = FatFraction(:,:,thorax_i);
%FF_abdomen = FatFraction(:,:,abdomen_i);
all_i = [thorax_i abdomen_i];
FF_organs = (organs3d(:,:,show_i)>0).*FatFraction(:,:,show_i);
FF_organs(FF_organs == 0) = NaN;
disp('--------------')
disp('Organs Fat Fraction: ')
disp('--------------')
X = ['Minimum           : ',num2str(min(FF_organs(:)))];
disp(X)
aveOrganFF = nanmean(FF_organs(:),1);
X = ['Average           : ',num2str(aveOrganFF)];
disp(X)
X = ['Maximum           : ',num2str(max(FF_organs(:)))];
disp(X)
if write_results
A = aveOrganFF;
xlRange = strcat('AA',num2str(row));
xlswrite(filename,A,sheet,xlRange)
end

disp('--------------')
disp('VAT volumes: ')
disp('--------------')
TATthorax = sum(TATvolume(thorax_i));
X = ['Thoracic           : ',num2str(TATthorax)];
disp(X)
VATabNoClip = sum(VATvolume(abdomen_no_clip_i));
X = ['Abdominal(no clip): ',num2str(VATabNoClip)];
disp(X)
VATabdomen = sum(VATvolume(show_i));
X = ['Heart max to FH:     ',num2str(VATabdomen)];
disp(X)
if write_results
A = [TATthorax VATabNoClip VATabdomen];
xlRange = strcat('AB',num2str(row));
xlswrite(filename,A,sheet,xlRange)
end

disp('--------------')
disp('IMAT volumes: ')
disp('--------------')
IMATthorax = sum(IMATvolume(thorax_i));
X = ['Thoracic           : ',num2str(IMATthorax)];
disp(X)
IMATabNoClip = sum(IMATvolume(abdomen_no_clip_i));
X = ['Abdominal(no clip): ',num2str(IMATabNoClip)];
disp(X)
IMATabdomen = sum(IMATvolume(abdomen_i));
X = ['Abdominal to FH:     ',num2str(IMATabdomen)];
disp(X)
if write_results
A = [IMATthorax IMATabNoClip IMATabdomen];
xlRange = strcat('AE',num2str(row));
xlswrite(filename,A,sheet,xlRange)
end

disp('---------------------')
disp('Slices of interest : ')
disp('---------------------')
EATsuperior_i = find(EATvolume > 0,1);
EATinferior_i =length(EATvolume) - find(fliplr(EATvolume) > 0,1) + 1;
EATsuperior_Z = upperTop + 1 - EATsuperior_i;
EATinferior_Z = upperTop + 1 - EATinferior_i;
% X = ['EAT superior i : ',num2str(EATsuperior_i)];
% disp(X)
% X = ['EAT inferior i : ',num2str(EATinferior_i)];
% disp(X)
X = ['EAT superior Z : ',num2str(EATsuperior_Z)];
disp(X)
X = ['EAT inferior Z : ',num2str(EATinferior_Z)];
disp(X)

PATsuperior_i = find(CATvolume > 0,1);
PATinferior_i =length(CATvolume) - find(fliplr(CATvolume) > 0,1) + 1;
PATsuperior_Z = upperTop + 1 - PATsuperior_i;
PATinferior_Z = upperTop + 1 - PATinferior_i;
% X = ['PAT superior i  : ',num2str(PATsuperior_i)];
% disp(X)
% X = ['PAT inferior i  : ',num2str(PATinferior_i)];
% disp(X)
X = ['PAT superior Z  : ',num2str(PATsuperior_Z)];
disp(X)
X = ['PAT inferior Z  : ',num2str(PATinferior_Z)];
disp(X)

PAATsuperior_i = find(PAATvolume > 0,1);
PAATinferior_i =length(PAATvolume) - find(fliplr(PAATvolume) > 0,1) + 1;
PAATsuperior_Z = upperTop + 1 - PAATsuperior_i;
PAATinferior_Z = upperTop + 1 - PAATinferior_i;
% X = ['PAAT superior i  : ',num2str(PAATsuperior_i)];
% disp(X)
% X = ['PAAT inferior i  : ',num2str(PAATinferior_i)];
% disp(X)
X = ['PAAT superior Z  : ',num2str(PAATsuperior_Z)];
disp(X)
X = ['PAAT inferior Z  : ',num2str(PAATinferior_Z)];
disp(X)

VATsuperior_i = find(VATvolume > 0,1);
VATsuperior_Z = upperTop + 1 - VATsuperior_i;
% X = ['VAT superior i : ',num2str(VATsuperior_i)];
% disp(X)
X = ['VAT superior Z : ',num2str(VATsuperior_Z)];
disp(X)

% X = ['SCAT superior i : ',num2str(thorax_i(1))];
% disp(X)
X = ['SCAT superior Z : ',num2str(upperTop)];
disp(X)

diaphram_Z = upperTop + 1 - thorax_i(end);
X = ['Diaphram (inferior of thorax) Z: ',num2str(diaphram_Z)];
disp(X)

X = ['Scan 1 superior Z  : ',num2str(upperTop)];
disp(X)
X = ['Scan 1 inferior Z  : ',num2str(upperBottom)];
disp(X)
X = ['Scan 2 superior Z  : ',num2str(lowerTop)];
disp(X) 
X = ['Scan 2 inferior Z  : ',num2str(lowerBottom)];
disp(X)

% 
% f1 = (3.0/10.5)*0.034
% disp('L1-L2 one slice estimator: ')
% (1/f1)*VATvolume(length(upperSlices))
% 
% f1 = (3.0/10.5)*0.041
% disp('L2-L3 one slice estimator: ')
% (1/f1)*0.5*(VATvolume(length(upperSlices) + Z_L1_L2 - Z_L2_L3) + VATvolume(end_i+3))
% 
% f1 = (3.0/10.5)*0.046
% disp('L3-L4 one slice estimator (best): ')
% (1/f1)*VATvolume(length(upperSlices) + Z_L1_L2 - Z_L3_L4) 
% 
% f1 = (3.0/10.5)*0.043
% disp('L4-L5 one slice estimator: ')
% (1/f1)*VATvolume(length(upperSlices) + Z_L1_L2 - Z_L4_L5) 
% 
% f1 = (3.0/10.5)*0.032
% disp('L5-S1 one slice estimator: ')
% (1/f1)*VATvolume(length(upperSlices) + Z_L1_L2 - Z_L5_S1b-3) 
% 
% f1 = (3.0/10.5)*0.04
% disp('UM one slice estimator: ')
% (1/f1)*VATvolume(length(upperSlices) + Z_L1_L2 - UM) 
% 
% f1 = (3.0/10.5)*0.014
% disp('FH one slice estimator: ')
% (1/f1)*VATvolume(end_i-2)
% 
% f5 = (3.0/10.5)*0.17
% disp('L1-L2 five slice estimator: ')
% (1/f5)*sum(VATvolume(length(upperSlices)-2:length(upperSlices)+2))
% 
% f5 = (3.0/10.5)*0.205
% disp('L2-L3 five slice estimator: ')
% (1/f5)*0.5*(sum(VATvolume(length(upperSlices) + Z_L1_L2 - Z_L2_L3-2:length(upperSlices) + Z_L1_L2 - Z_L2_L3+2)) ...
%           +sum(+ VATvolume(end_i+1:end_i+5)))
%       
% f5 = (3.0/10.5)*0.226
% disp('L3-L4 five slice estimator (best): ')
% (1/f5)*sum(VATvolume(length(upperSlices) + Z_L1_L2 - Z_L3_L4-2:length(upperSlices) + Z_L1_L2 - Z_L3_L4+2))
% 
% f5 = (3.0/10.5)*0.209
% disp('L4-L5 five slice estimator: ')
% (1/f5)*sum(VATvolume(length(upperSlices) + Z_L1_L2 - Z_L4_L5-2:length(upperSlices) + Z_L1_L2 - Z_L4_L5+2))
% 
% f5 = (3.0/10.5)*0.162
% disp('L5-S1 five slice estimator: ')
% (1/f5)*sum(VATvolume(length(upperSlices) + Z_L1_L2 - Z_L5_S1b-3-2:length(upperSlices) + Z_L1_L2 - Z_L5_S1b-3+2))
% 
% f5 = (3.0/10.5)*0.198
% disp('UM five slice estimator: ')
% (1/f5)*sum(VATvolume(length(upperSlices) + Z_L1_L2 - UM-2:length(upperSlices) + Z_L1_L2 - UM+2))
% 
% f5 = (3.0/10.5)*0.104
% disp('FH five slice estimator: ')
% (1/f5)*sum(VATvolume(end_i-4:end_i))
% 

if showSlice
% render a sagittal slice   
    slice2show = uint8(zeros(size(aorta3d,1),size(aorta3d,2),3));
    % show aorta as magenta
    slice2show(:,:,1) = slice2show(:,:,1) + aorta3d(:,:,sliceNbr);
    slice2show(:,:,3) = slice2show(:,:,3) + aorta3d(:,:,sliceNbr);    
    % show PAAT as cyan
    slice2show(:,:,2) = slice2show(:,:,2) + PAAT3d(:,:,sliceNbr);
    slice2show(:,:,3) = slice2show(:,:,3) + PAAT3d(:,:,sliceNbr);   
    % show heart as red
    slice2show(:,:,1) = slice2show(:,:,1) + heart3d(:,:,sliceNbr);    
    % show EAT as silver
    slice2show(:,:,1) = slice2show(:,:,1) + uint8(EAT3d(:,:,sliceNbr));    
    slice2show(:,:,2) = slice2show(:,:,2) + uint8(EAT3d(:,:,sliceNbr));
    slice2show(:,:,3) = slice2show(:,:,3) + uint8(EAT3d(:,:,sliceNbr));  
    % show CAT (PAT) as sea green   
    slice2show(:,:,2) = slice2show(:,:,2) + uint8(CAT3d(:,:,sliceNbr));
    slice2show(:,:,3) = slice2show(:,:,3) + uint8(0.667*CAT3d(:,:,sliceNbr));
%    % show voids as dark blue  
%    slice2show(:,:,3) = slice2show(:,:,3) + uint8(0.5*voids3d(:,:,sliceNbr));  
%    % show lungs blue  
%    slice2show(:,:,3) = slice2show(:,:,3) + uint8(lung3d(:,:,sliceNbr)); 
    % show TAT as dark green  
    slice2show(:,:,2) = slice2show(:,:,2) + uint8(0.667*TAT3d(:,:,sliceNbr));         
    % show VAT as green  
    slice2show(:,:,2) = slice2show(:,:,2) + uint8(VAT3d(:,:,sliceNbr)); 
    % show organs as dark red  
    slice2show(:,:,1) = slice2show(:,:,1) + uint8(0.5*organs3d(:,:,sliceNbr)); 
    % show IMAT as beige 
    slice2show(:,:,1) = slice2show(:,:,1) + uint8(IMAT3d(:,:,sliceNbr)); 
    slice2show(:,:,2) = slice2show(:,:,2) + uint8(0.75*IMAT3d(:,:,sliceNbr)); 
    slice2show(:,:,3) = slice2show(:,:,3) + uint8(0.5*IMAT3d(:,:,sliceNbr)); 
    % show SCAT as yellow
    slice2show(:,:,1) = slice2show(:,:,1) + uint8(SCAT3d(:,:,sliceNbr));    
    slice2show(:,:,2) = slice2show(:,:,2) + uint8(SCAT3d(:,:,sliceNbr)); 
    slice2show(:,:,3) = slice2show(:,:,3) + uint8(0.333*SCAT3d(:,:,sliceNbr)); 
    
    figure, imshow(slice2show);
end

if strcmp(render,'VAT')
% 3D render VAT
VAT3d = squeeze(VAT3d(:,:,1:end_i));
VAT3d = flip(VAT3d,3);
figure,
hv = vol3d('cdata',VAT3d ,'texture','3D');
view(3);  
axis([50 275 50 200 0 inf]);  daspect([1 1 1.40625/3.0])
colormap(autumn(256));
alphamap('rampup');
% alphamap(.06 .* alphamap);
% %alphamap([0 linspace(0.1, 0, 255)]);
end

if strcmp(render,'organs')
% 3D render organs
organs3d = squeeze(organs3d(:,:,1:end_i));
organs3d = flip(organs3d,3);
figure,
ho = vol3d('cdata',organs3d ,'texture','3D');
view(3);  
axis([75 250 75 175 0 inf]);  daspect([1 1 1.40625/3.0])
colormap(copper(256));
%alphamap('rampup');
alphamap(.06 .* alphamap);
alphamap([0 linspace(0.1, 0, 255)]);
end

if strcmp(render,'lung')
% 3D render lung
lung3d = squeeze(lung3d);
lung3d = flip(lung3d,3);
figure,
hl = vol3d('cdata',lung3d ,'texture','3D');
view(3);  
axis([50 275 50 200 0 inf]);  daspect([1 1 1.40625/3.0])
colormap(winter(256));
alphamap('rampup');
%alphamap(.36 .* alphamap);
%alphamap([0 linspace(0.1, 0, 255)]);
end

if strcmp(render,'heart')
% 3D render heart
heart3d = squeeze(heart3d);
heart3d = flip(heart3d,3);
figure,
hh = vol3d('cdata',heart3d ,'texture','3D');
view(3);  
axis([75 200 50 175 0 inf]);  daspect([1 1 1.40625/3.0])
colormap(copper(256));
alphamap('rampup');
% alphamap(.06 .* alphamap);
% %alphamap([0 linspace(0.1, 0, 255)]);
end

if strcmp(render,'CAT') || strcmp(render,'PAT')
% 3D render CAT
CAT3d = squeeze(CAT3d);
CAT3d = flip(CAT3d,3);
figure,
hc = vol3d('cdata',CAT3d ,'texture','3D');
view(3);  
axis([75 200 50 175 0 inf]);   daspect([1 1 1.40625/3.0])
colormap(autumn(256));
alphamap('rampup');
% alphamap(.06 .* alphamap);
% alphamap([0 linspace(0.1, 0, 255)]);
end

if strcmp(render,'EAT')
% 3D render EAT
EAT3d = squeeze(EAT3d);
EAT3d = flip(EAT3d,3);
figure,
hc = vol3d('cdata',EAT3d ,'texture','3D');
view(3);  
axis([75 200 50 175 0 inf]);   daspect([1 1 1.40625/3.0])
colormap(autumn(256));
alphamap('rampup');
% alphamap(.06 .* alphamap);
% alphamap([0 linspace(0.1, 0, 255)]);
end

if strcmp(render,'aorta')
% 3D render aorta
aorta3d = squeeze(aorta3d);
aorta3d = flip(aorta3d,3);
figure,
ha = vol3d('cdata',aorta3d ,'texture','3D');
view(3);  
axis([125 175 125 175 0 inf]);   daspect([1 1 1.40625/3.0])
colormap(cool(256));
alphamap('rampup');
% alphamap(.06 .* alphamap);
end

if strcmp(render,'PAAT')
% 3D render PAAT
PAAT3d = squeeze(PAAT3d);
PAAT3d = flip(PAAT3d,3);
figure,
hp = vol3d('cdata',PAAT3d ,'texture','3D');
view(3);  
axis([125 175 125 175 0 inf]);   daspect([1 1 1.40625/3.0])
colormap(autumn(256));
alphamap('rampup');
% alphamap(.06 .* alphamap);
end



