function [br,T] = add_roi(br,roi_size,roi_start,T,y)

z_size = roi_size(1)+1;
x_size = roi_size(2)+1;
num_slices = roi_size(3)+1;
z_start = roi_start(1);
x_start = roi_start(2);
slice_start = roi_start(3);
[x, yy, z] = meshgrid(1:z_size, 1:x_size, 1:num_slices);
r = sqrt((x-z_size/2).^2 + (yy-x_size/2).^2 + (z-num_slices/2).^2);
r = r(1:x_size-1,1:z_size-1,1:num_slices-1);
r = max(r(:)) - r;
r = r - min(r(:));
br(x_start:(x_start+x_size-2),z_start:(z_start+z_size-2),...
    slice_start:(slice_start+num_slices-2)) = r/4;
for i = slice_start:(slice_start + num_slices - 2)
    for j = z_start:(z_start+z_size-2)
        for k = x_start:(x_start+x_size-2)
            T(j,k,i,:) = squeeze(T(j,k,i,:)) + br(j,k,i)*y;
        end
    end
end

end