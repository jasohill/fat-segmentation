% Author: Jason E. Hill
% updated 3 MAY 2017
function outVol = translate_lower(inVol, Nupper, shift)
outVol = inVol;

temp = inVol;
temp(Nupper) = temp(Nupper+1);

x = 1:length(temp);
v = double(squeeze(temp(:)));
xq = x + shift; 
vq = interp1(x,v,xq);

outVol(Nupper+1:end) = vq(Nupper+1:end);

end

