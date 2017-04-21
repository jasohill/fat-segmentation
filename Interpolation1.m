% 

%% Interpolation Grid

W = 33; % Width of intravertaal mask
H = 26; % Height of mask 

x = 159; % Middle of disk x = [143 X 177]; % Dimensions of intevertabral disk
y = 126; % Middle of disk y = [112 Y 138]; % Dimensions of intervertabral disk

DW = 320; % Slice Size
DH = 260; % Slice Size

%deltaZ = [Z-1 Z Z+1]
%I = [ ];
%p = polyfit(deltaZ, I, 2);
% z_peak = -b/(2a);

[X Y] = meshgrid(1:DW,1:DH); % Image size 
EMask = ((X - x)/W).^2 + ((Y - y)/H).^2 <= 1;

% imagesc(EMask);