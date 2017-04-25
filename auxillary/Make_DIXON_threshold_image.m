function Bdata = Make_DIXON_threshold_image(Fdata, Wdata)
% Creates a combined threshold image from fat and water saturated DIXON
% data.

dims = size(Fdata);

maxData = min([min((max(max(Fdata)))),min((max(max(Wdata))))]);

Fdata(Fdata>maxData) = maxData;
Wdata(Wdata>maxData) = maxData;
maxData1 = double(maxData);

for slice = 1:dims(3)
    Bdata(:,:,slice) = double(fliplr(rot90(Fdata(:,:,slice)',2)))/maxData1 ...
                  + double(fliplr(rot90(Wdata(:,:,slice)',2)))/maxData1;          
end
maxData2 = min([maxData1,min(max(max(Bdata)))*maxData1]);
Bdata = Bdata*(maxData1/maxData2);
Bdata(Bdata>1.0) = 1.0;
Bdata(Bdata<0.1) = 0.0;

end

