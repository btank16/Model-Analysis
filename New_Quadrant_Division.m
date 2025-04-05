close all
clear all

% Select files of interest
[file,path] = uigetfile('*.nii.gz','Multiselect','on');
if ~iscell(file)
    file = {file};
end

    clearvars -except file_count file path DataTable
    %% Determine file name
    niftifile = char(file);
    filename_tmp = strsplit(niftifile,'.');
    filename = filename_tmp{1};
    lumen_stl = strcat(filename,'_Lumen.stl');% Name for lumen stl file 
    wall_stl  = strcat(filename,'_Wall.stl'); % Name of wall stl file 
    thrombus_stl = strcat(filename, '_lumen_outer.stl'); % Name of thrombus stl file

    %% AQUIRE ULTRASOUND INFO NEEDED
    %Important: This code relies on isotropic voxels
    info = niftiinfo([path niftifile]);  %retrieve info from ultraound NifTi file
    PixelDimensions = info.PixelDimensions; % Gives a matrix of voxel size: should be isotropic
    voxel_size = PixelDimensions(:,1);
    %% READ STL FILE 
    [vert_lumen_d,faces_lumen_d,normals_lumen_d,name_lumen_d] = stlRead([path lumen_stl]);
    [vert_wall_d,faces_wall_d,normals_wall_d,name_wall_d] = stlRead([path wall_stl]);
    [vert_thrombus_d,faces_thrombus_d,normals_thrombus_d,name_thrombus_d] = stlRead([path thrombus_stl]);

    %% CONVERT STL FILE TO A SOLID MESH

    %Create grids 
    x_ultrasound_grid = linspace(0, (info.ImageSize(2)-1)*voxel_size, info.ImageSize(2));
    y_ultrasound_grid = linspace(0, (info.ImageSize(1)-1)*voxel_size, info.ImageSize(1));
    z_ultrasound_grid = linspace(0, (info.ImageSize(3)-1)*voxel_size, info.ImageSize(3));

    % Create a lumen structure array containing 2 field variables: (1)the vertices and (2) faces 
    lumen_d_FV(1).vertices = vert_lumen_d;
    lumen_d_FV(1).faces = faces_lumen_d;

    % Create a wall structure array containing 2 field variables: (1)the vertices and (2) faces 
    wall_d_FV(1).vertices = vert_wall_d;
    wall_d_FV(1).faces = faces_wall_d;

    % Create a thrombus structure array containing 2 field variables: (1)the vertices and (2) faces
    thrombus_d_FV(1).vertices = vert_thrombus_d;
    thrombus_d_FV(1).faces = faces_thrombus_d;
    
    % Convert STL to a Solid Mesh
    [lumen_d_mesh] = VOXELISE(y_ultrasound_grid,x_ultrasound_grid,z_ultrasound_grid,lumen_d_FV);
    [wall_d_mesh] = VOXELISE(y_ultrasound_grid,x_ultrasound_grid,z_ultrasound_grid,wall_d_FV);
    [thrombus_d_mesh] = VOXELISE(y_ultrasound_grid,x_ultrasound_grid,z_ultrasound_grid,thrombus_d_FV);

    % Fix Orientation for Flipping of X and Y between Matrix and DataPoints
    lumen_d_mesh = permute(lumen_d_mesh,[2 1 3]);
    wall_d_mesh = permute(wall_d_mesh,[2 1 3]);
    thrombus_d_mesh = permute(thrombus_d_mesh, [2 1 3]);
    
    
z_length = length(z_ultrasound_grid);
centerline = zeros(z_length,3);


for a = 1:z_length
    if sum(thrombus_d_mesh(:,:,a),'all')~=0
    z_slice = regionprops(thrombus_d_mesh(:,:,a),'Centroid');
    if isequal(size(z_slice),[1,1])
    x_center = z_slice.Centroid(:,1);
    y_center = z_slice.Centroid(:,2);
    centerline(a,1) = x_center;
    centerline(a,2) = y_center;
    end
    end
end
    
centerline(:,3) = 1:z_length;    

 y_grid_max = length(y_ultrasound_grid);
 x_grid_max = length(x_ultrasound_grid);
 Q1_mesh = zeros(x_grid_max,y_grid_max,length(centerline));
 Q2_mesh = zeros(x_grid_max,y_grid_max,length(centerline));
 Q3_mesh = zeros(x_grid_max,y_grid_max,length(centerline));
 Q4_mesh = zeros(x_grid_max,y_grid_max,length(centerline));


for b = 1:length(centerline)
    x_limit = round(centerline(b,1));
    y_limit = round(centerline(b,2));
    c = centerline(b,3);
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
    if x_limit ~= 0
    Q4_mesh(ylim_upper:x_grid_max,xlim_upper:y_grid_max,b) = thrombus_d_mesh(ylim_upper:x_grid_max,xlim_upper:y_grid_max,c);
    Q3_mesh(ylim_upper:x_grid_max,1:xlim_lower,b) = thrombus_d_mesh(ylim_upper:x_grid_max,1:xlim_lower,b);
    Q2_mesh(1:ylim_lower,xlim_upper:y_grid_max,b) = thrombus_d_mesh(1:ylim_lower,xlim_upper:y_grid_max,b);
    Q1_mesh(1:ylim_lower,1:xlim_lower,b) = thrombus_d_mesh(1:ylim_lower,1:xlim_lower,c);
    end
end

Q1_area = zeros(z_length,2);
Q2_area = zeros(z_length,2);
Q3_area = zeros(z_length,2);
Q4_area = zeros(z_length,2);

for d = 1:z_length
    % Area calculation for Q1
    Q1_pixel_area = sum(Q1_mesh(:,:,d),'all');
    Q1_area(d,1) = Q1_pixel_area;
    Q1_area(d,2) = (Q1_pixel_area*voxel_size);
    % Area calculation for Q2
    Q2_pixel_area = sum(Q2_mesh(:,:,d),'all');
    Q2_area(d,1) = Q2_pixel_area;
    Q2_area(d,2) = (Q2_pixel_area*voxel_size);
    %Area calculation for Q3
    Q3_pixel_area = sum(Q3_mesh(:,:,d),'all');
    Q3_area(d,1) = Q3_pixel_area;
    Q3_area(d,2) = (Q3_pixel_area*voxel_size);
    % Area calculation for Q4
    Q4_pixel_area = sum(Q4_mesh(:,:,d),'all');
    Q4_area(d,1) = Q4_pixel_area;
    Q4_area(d,2) = (Q4_pixel_area*voxel_size);
end
    
    
    
    
    
    
    
    
