function out = map2original(input)
% input: 374 x k data.
% mask. Combined_Parcels.nii

maskobj = fmri_data(which('Combined_Parcels.nii'));
out = maskobj;
out.dat = zeros(size(maskobj.dat, 1), size(input, 2));
for i = 1:size(input, 2)
    for j = 1:max(maskobj.dat)
        out.dat(maskobj.dat == j, i) = input(j, i);
    end
end


end