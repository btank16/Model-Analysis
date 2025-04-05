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
    
%Call models into the Mesh_Quadrants function
[lumen_centerline,Q1_lumen,Q2_lumen,Q3_lumen,Q4_lumen,Q1_lumen_area,Q2_lumen_area,Q3_lumen_area,Q4_lumen_area] = Mesh_Quadrants(lumen_d_mesh,voxel_size,x_ultrasound_grid,y_ultrasound_grid,z_ultrasound_grid);
[thrombus_centerline,Q1_thrombus,Q2_thrombus,Q3_thrombus,Q4_thrombus,Q1_thrombus_area,Q2_thrombus_area,Q3_thrombus_area,Q4_thrombus_area] = Mesh_Quadrants(thrombus_d_mesh,voxel_size,x_ultrasound_grid,y_ultrasound_grid,z_ultrasound_grid);

%Subtracts the lumen area from the lumen_outer area to get thrombus area
Q1_thrombusA = Q1_thrombus_area(:,2)-Q1_lumen_area(:,2);  
Q2_thrombusA = Q2_thrombus_area(:,2)-Q2_lumen_area(:,2);  
Q3_thrombusA = Q3_thrombus_area(:,2)-Q3_lumen_area(:,2);  
Q4_thrombusA = Q4_thrombus_area(:,2)-Q4_lumen_area(:,2);  

%Removes "negative" thrombus areas
for a = 1:length(z_ultrasound_grid)
    if Q1_thrombusA(a)<0
        Q1_thrombusA(a) = 0;
    end
    if Q2_thrombusA(a)<0
        Q2_thrombusA(a) = 0;
    end
    if Q3_thrombusA(a)<0
        Q3_thrombusA(a) = 0;
    end
    if Q4_thrombusA(a)<0
        Q4_thrombusA(a) = 0;
    end
end
    
%Plots the four quadrants along the length of the vessel
plot(Q1_thrombus_area(:,3),Q1_thrombusA,'.r');
hold on
plot(Q2_thrombus_area(:,3),Q2_thrombusA,'.g');
plot(Q3_thrombus_area(:,3),Q3_thrombusA,'.b');
plot(Q4_thrombus_area(:,3),Q4_thrombusA,'.k');
hold off

%Calculates thrombus volume
Q1_burden = (sum(Q1_thrombusA)*voxel_size);
Q2_burden = (sum(Q2_thrombusA)*voxel_size);
Q3_burden = (sum(Q3_thrombusA)*voxel_size);
Q4_burden = (sum(Q4_thrombusA)*voxel_size);
total_burden = Q1_burden+Q2_burden+Q3_burden+Q4_burden;

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    