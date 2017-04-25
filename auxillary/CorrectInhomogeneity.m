function [J,Icorrect] = CorrectInhomogeneity(I,radius,varargin)
% construct rough estimate of the signal inhomogeneity via morphological closing

factor = inf;
% handle variable inputs
if nargin == 4
    factor = varargin{2};
    filter = varargin{1};
elseif nargin == 3
    filter = varargin{1};
else
    filter = [];
end

% estimate inhomogeneity field
Iinhomogeneity = imclose(I,strel('disk',radius));

% imshow(Iinhomogeneity,[])

% construct the correction field by the estimated signal inhomogeneity
Icorrection = (mean(mean(Iinhomogeneity)))./double(Iinhomogeneity);
Icorrection(Icorrection>factor) = 1./Icorrection(Icorrection>factor);
%figure, imshow(Icorrection,[])
Icorrect = double(I).*Icorrection;
if ~isempty(filter)
    if strcmp(filter,'medfilt')
        Icorrect = medfilt2(Icorrect,[3,3]);
    else
        warning('unsupported option for filter provided')
    end
end
% rescale image to the maximum intensity
Imax = max(max(Icorrect));
J = uint8((Icorrect/Imax)*255);
end

