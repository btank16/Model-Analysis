function [path_length,vessel_length,tortuosity] = Tortuosity_Function(mesh,voxel_size,z_ultrasound_grid) 
    
z_length = length(z_ultrasound_grid);
centerline = zeros(z_length,2);

%Finds the centerline of the stl model
for a = 1:z_length
    if sum(mesh(:,:,a),'all')~=0
    z_slice = regionprops(mesh(:,:,a),'Centroid');
    if isequal(size(z_slice),[1,1])
    x_center = z_slice.Centroid(:,1);
    y_center = z_slice.Centroid(:,2);
    centerline(a,1) = x_center;
    centerline(a,2) = y_center;
    end
    end
end
    
deleteZrows = centerline(:,1)<=0;
centerline(deleteZrows,:) = [];
    
point_distance = zeros(length(centerline),1);

%Finds the length between each point in the centerline 
for b = 1:(length(centerline)-1)
    %Finds the distance apart in the xy-plane
    x_length = abs(centerline(b,1)-centerline(b+1,1));
    y_length = abs(centerline(b,2)-centerline(b+1,2));
    xy_length = sqrt((x_length^2)+(y_length^2))*voxel_size;
    %finds distance between points
    xyz_length = sqrt((xy_length^2)+(voxel_size^2));
    point_distance(b,1) = xyz_length;
end
  
path_length = sum(point_distance);
vessel_length = (length(centerline)*voxel_size);
tortuosity = path_length/vessel_length;
    
end    
    
    
    
    
    
    
    
    
    
    