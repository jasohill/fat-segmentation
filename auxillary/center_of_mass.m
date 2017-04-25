function CenterOfMass = center_of_mass(data)
% Find the center of mass of a 3d or 2d data array.

   size_data = size(data);
   dims = length(size_data);
   
   mass = double(sum(data(:)));
   CenterOfMass = zeros(dims,1);

   for dim = 1:dims
       size_dim = ones(1,dims);
       size_dim(dim) = size_data(dim);
       replicate = size_data;
       replicate(dim) = 1;
       indices = repmat(reshape(1:size_data(dim),size_dim),replicate);
       CenterOfMass(dim) = sum(double(indices(:)).*double(data(:)))./mass;
   end

end

