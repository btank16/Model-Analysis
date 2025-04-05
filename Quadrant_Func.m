function [Q1_mesh,Q2_mesh,Q3_mesh,Q4_mesh,Q1_area,Q2_area,Q3_area,Q4_area] = Quadrant_Func(mesh,centerline,voxel_size,x_ultrasound_grid,y_ultrasound_grid,z_ultrasound_grid)

z_length = length(z_ultrasound_grid);
 y_grid_max = length(y_ultrasound_grid);
 x_grid_max = length(x_ultrasound_grid);
 Q1_mesh = zeros(x_grid_max,y_grid_max,length(centerline));
 Q2_mesh = zeros(x_grid_max,y_grid_max,length(centerline));
 Q3_mesh = zeros(x_grid_max,y_grid_max,length(centerline));
 Q4_mesh = zeros(x_grid_max,y_grid_max,length(centerline));

%Divides the image in four quadrants at each z-slice
for b = 1:length(centerline)
    x_limit = round(centerline(b,1));
    y_limit = round(centerline(b,2));
    c = centerline(b,3);
    % Sets upper and lower bounds to avoid overlap
    if x_limit > centerline(b,1)
        xlim_upper = x_limit;
        xlim_lower = x_limit - 1;
    else
        xlim_upper = x_limit + 1;
        xlim_lower = x_limit;
    end
    if y_limit > centerline(b,2)
        ylim_upper = y_limit;
        ylim_lower = y_limit - 1;
    else
        ylim_upper = y_limit + 1;
        ylim_lower = y_limit;
    end
    %Places each region of interest in each quadrant
    if x_limit ~= 0
    Q4_mesh(ylim_upper:x_grid_max,xlim_upper:y_grid_max,b) = mesh(ylim_upper:x_grid_max,xlim_upper:y_grid_max,c);
    Q3_mesh(ylim_upper:x_grid_max,1:xlim_lower,b) = mesh(ylim_upper:x_grid_max,1:xlim_lower,b);
    Q2_mesh(1:ylim_lower,xlim_upper:y_grid_max,b) = mesh(1:ylim_lower,xlim_upper:y_grid_max,b);
    Q1_mesh(1:ylim_lower,1:xlim_lower,b) = mesh(1:ylim_lower,1:xlim_lower,c);
    end
end

Q1_area = zeros(z_length,3);
Q2_area = zeros(z_length,3);
Q3_area = zeros(z_length,3);
Q4_area = zeros(z_length,3);

%Calculates the area at each z-slice
for d = 1:z_length
    % Area calculation for Q1
    Q1_pixel_area = sum(Q1_mesh(:,:,d),'all');
    Q1_area(d,1) = Q1_pixel_area;
    Q1_area(d,2) = (Q1_pixel_area*voxel_size^2);
    % Area calculation for Q2
    Q2_pixel_area = sum(Q2_mesh(:,:,d),'all');
    Q2_area(d,1) = Q2_pixel_area;
    Q2_area(d,2) = (Q2_pixel_area*voxel_size^2);
    %Area calculation for Q3
    Q3_pixel_area = sum(Q3_mesh(:,:,d),'all');
    Q3_area(d,1) = Q3_pixel_area;
    Q3_area(d,2) = (Q3_pixel_area*voxel_size^2);
    % Area calculation for Q4
    Q4_pixel_area = sum(Q4_mesh(:,:,d),'all');
    Q4_area(d,1) = Q4_pixel_area;
    Q4_area(d,2) = (Q4_pixel_area*voxel_size^2);
end

Q1_area(:,3) = z_ultrasound_grid.';
Q2_area(:,3) = z_ultrasound_grid.';
Q3_area(:,3) = z_ultrasound_grid.';
Q4_area(:,3) = z_ultrasound_grid.';


end