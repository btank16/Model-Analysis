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
    thrombus_mesh = thrombus_d_mesh-lumen_d_mesh;
    
    
%Call models into the Mesh_Quadrants function
% [lumen_centerline,Q1_lumen,Q2_lumen,Q3_lumen,Q4_lumen,Q1_lumen_area,Q2_lumen_area,Q3_lumen_area,Q4_lumen_area] = Mesh_Quadrants(lumen_d_mesh,voxel_size,x_ultrasound_grid,y_ultrasound_grid,z_ultrasound_grid);
% [thrombus_centerline,Q1_thrombus,Q2_thrombus,Q3_thrombus,Q4_thrombus,Q1_thrombus_area,Q2_thrombus_area,Q3_thrombus_area,Q4_thrombus_area] = Mesh_Quadrants(thrombus_d_mesh,voxel_size,x_ultrasound_grid,y_ultrasound_grid,z_ultrasound_grid);
%[centerline,Q1_thrombusSET,Q2_thrombusSET,Q3_thrombusSET,Q4_thrombusSET,Q1_thrombusA,Q2_thrombusA,Q3_thrombusA,Q4_thrombusA] = Mesh_Quadrants(thrombus_mesh,voxel_size,x_ultrasound_grid,y_ultrasound_grid,z_ultrasound_grid);

%Find Centerline of each model
[lumen_center] = Centroid_Func(lumen_d_mesh,z_ultrasound_grid);
[thrombus_center] = Centroid_Func(thrombus_d_mesh,z_ultrasound_grid);

%Splits mesh into quadrants using the lumen centerline
[Q1_lumen,Q2_lumen,Q3_lumen,Q4_lumen,Q1_lumen_area,Q2_lumen_area,Q3_lumen_area,Q4_lumen_area] = Quadrant_Func(lumen_d_mesh,lumen_center,voxel_size,x_ultrasound_grid,y_ultrasound_grid,z_ultrasound_grid);
[Q1_thrombus,Q2_thrombus,Q3_thrombus,Q4_thrombus,Q1_thrombus_area,Q2_thrombus_area,Q3_thrombus_area,Q4_thrombus_area] = Quadrant_Func(thrombus_d_mesh,lumen_center,voxel_size,x_ultrasound_grid,y_ultrasound_grid,z_ultrasound_grid);

%Subtracts the lumen area from the lumen_outer area to get thrombus area
Q1_thrombusA = Q1_thrombus_area(:,2)-Q1_lumen_area(:,2);  
Q2_thrombusA = Q2_thrombus_area(:,2)-Q2_lumen_area(:,2);  
Q3_thrombusA = Q3_thrombus_area(:,2)-Q3_lumen_area(:,2);  
Q4_thrombusA = Q4_thrombus_area(:,2)-Q4_lumen_area(:,2);  

%Flips the mesh to start in proximal position
Q1_thrombusA = flip(Q1_thrombusA);
Q2_thrombusA = flip(Q2_thrombusA);
Q3_thrombusA = flip(Q3_thrombusA);
Q4_thrombusA = flip(Q4_thrombusA);

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

total_thrombusA = Q1_thrombusA+Q2_thrombusA+Q3_thrombusA+Q4_thrombusA;


%Plots the four quadrants along the length of the vessel
% plot(Q1_thrombus_area(:,3),Q1_thrombusA,'.r');
% hold on
% plot(Q2_thrombus_area(:,3),Q2_thrombusA,'.g');
% plot(Q3_thrombus_area(:,3),Q3_thrombusA,'.b');
% plot(Q4_thrombus_area(:,3),Q4_thrombusA,'.k');
% % plot(Q4_thrombus_area(:,3),total_thrombusA,'.m');
% xlabel('Vessel Length (mm)')
% ylabel('Thrombus Cross-Sectional Area (mm^2)')
% title('Thrombus Area per Quadrant Along the Length of the Vessel')
% hold off
% legend('Dorsal/Right','Dorsal/Left','Ventral/Right','Ventral/Left')

%Calculates thrombus volume
Q1_burden = (sum(Q1_thrombusA)*voxel_size);
Q2_burden = (sum(Q2_thrombusA)*voxel_size);
Q3_burden = (sum(Q3_thrombusA)*voxel_size);
Q4_burden = (sum(Q4_thrombusA)*voxel_size);
total_burden = Q1_burden+Q2_burden+Q3_burden+Q4_burden;

%Call Tortuosity Function
[lumen_path,lumen_length,lumen_tortuosity] = Tortuosity_Function(lumen_d_mesh,voxel_size,z_ultrasound_grid);
[thrombus_path,thrombus_length,thrombus_tortuosity] = Tortuosity_Function(thrombus_d_mesh,voxel_size,z_ultrasound_grid);

%Calculate ration of lumen_tortuosity/thrombus_tortuosity  
tortuosity_ratio = lumen_tortuosity/thrombus_tortuosity;
    
%Creates a table to organize the output variables    
% datatable = [{'Q1 Burden'} Q1_burden; {'Q2 Burden'} Q2_burden;{'Q3 Burden'} Q3_burden;{'Q4 Burden'} Q4_burden;{'Total Burden'} total_burden;{'Lumen Tortuosity'} lumen_tortuosity;{'Thrombus Tortuosity'} thrombus_tortuosity;{'Tortuosity Ratio'} tortuosity_ratio]
    
%Finding percent area per quadrant along the vessel length
z_length = length(z_ultrasound_grid);
percent_area = zeros(z_length,4);
for b = 1:z_length
    if total_thrombusA(b,1)~=0
        percent_area(b,1) = (Q1_thrombusA(b,1)/total_thrombusA(b,1))*100;
        percent_area(b,2) = (Q2_thrombusA(b,1)/total_thrombusA(b,1))*100;
        percent_area(b,3) = (Q3_thrombusA(b,1)/total_thrombusA(b,1))*100;
        percent_area(b,4) = (Q4_thrombusA(b,1)/total_thrombusA(b,1))*100;
    end
end

%Finds the maximum lumen and thrombus with position
vessel_area = Q1_thrombus_area(:,2)+Q2_thrombus_area(:,2)+Q3_thrombus_area(:,2)+Q4_thrombus_area(:,2);
vessel_area = flip(vessel_area);
[max_area,max_vessel] = max(vessel_area);
[max_areaT,max_thrombus] = max(total_thrombusA);
vessel_distance = max_vessel*voxel_size;
thrombus_dist = max_thrombus*voxel_size;
vessel_length = z_length*voxel_size;
 datatable2 = [{'Max Area'} max_area;{'Vessel Distance'} vessel_distance;{'Max Thrombus'} max_areaT;{'Thrombus Distance'} thrombus_dist;{'Vessel Length'} vessel_length]

%Removes any values of 100% from quadrant matrix
for c = 1:length(z_ultrasound_grid)
    if percent_area(c,1)==100
        percent_area(c,1) = 0;
    end
    if percent_area(c,2)==100
        percent_area(c,2) = 0;
    end
    if percent_area(c,3)==100
        percent_area(c,3) = 0;
    end
    if percent_area(c,4)==100
        percent_area(c,4) = 0;
    end
end

%Finds the maximum percentage for each quadrant and its position
[maxQ1,Q1a] = max(percent_area(:,1));
[maxQ2,Q2a] = max(percent_area(:,2));
[maxQ3,Q3a] = max(percent_area(:,3));
[maxQ4,Q4a] = max(percent_area(:,4));
Q1_dist = Q1a*voxel_size;
Q2_dist = Q2a*voxel_size;
Q3_dist = Q3a*voxel_size;
Q4_dist = Q4a*voxel_size;
% datatable3 = [{'Q1 %'} maxQ1;{'Q1 Dist'} Q1_dist;{'Q2 %'} maxQ2;{'Q2 Dist'} Q2_dist;{'Q3 %'} maxQ3;{'Q3 Dist'} Q3_dist;{'Q4 %'} maxQ4;{'Q4 Dist'} Q4_dist]



%Plot percent area over the length of the vessel
% plot(Q1_thrombus_area(:,3),percent_area(:,1),'.r');
% hold on
% plot(Q2_thrombus_area(:,3),percent_area(:,2),'.g');
% plot(Q3_thrombus_area(:,3),percent_area(:,3),'.b');
% plot(Q4_thrombus_area(:,3),percent_area(:,4),'.k');
% hold off
% legend('Dorsal/Right','Dorsal/Left','Ventral/Right','Ventral/Left')

lumen_centerX = lumen_center(:,1);
lumen_centerY = lumen_center(:,2);
thrombus_centerX = thrombus_center(:,1);
thrombus_centerY = thrombus_center(:,2);

centroid_diff = zeros(length(lumen_center),2);
%Find distance and angle between 
for d = 1:length(lumen_center)
    x_len = -(thrombus_centerX(d)-lumen_centerX(d))*voxel_size;
    y_len = -(thrombus_centerY(d)-lumen_centerY(d))*voxel_size;
    cent_diff = sqrt((abs(x_len)^2)+(abs(y_len)^2));
    if cent_diff>4
        cent_diff = 0;
    end
    centroid_diff(d,1) = cent_diff;
    if cent_diff ~=0
        if x_len>0 && y_len>0
            cent_ang = asin(y_len/cent_diff);
            centroid_diff(d,2) = cent_ang;
        elseif x_len==0 && y_len>0
            centroid_diff(d,2) = pi/2;
        elseif x_len==0 && y_len<0
            centroid_diff(d,2) = (3*pi)/2;
        elseif x_len>0 && y_len==0
            centroid_diff(d,2) = 0;
        elseif x_len<0 && y_len==0
            centroid_diff(d,2) = pi;
        elseif x_len>0 && y_len<0
            cent_ang = (2*pi)-acos(x_len/cent_diff);
            centroid_diff(d,2) = cent_ang;
        elseif x_len<0 && y_len>0
            cent_ang = pi-asin(y_len/cent_diff);
            centroid_diff(d,2) = cent_ang;
        elseif x_len<0 && y_len<0
            cent_ang = atan(y_len/x_len)+pi;
            centroid_diff(d,2) = cent_ang;
        end
    end
end

%finds max difference between centroids and the angle at that point
[max_Centdiff,Cent_dist] = max(centroid_diff(:,1));
% max_ang = centroid_diff(Cent_dist,2);
total_AREA = flip(total_thrombusA);
[max_T,max_AREA] = max(total_AREA);

% imshow(thrombus_mesh(:,:,max_AREA))
% hold on
% plot(lumen_centerX(max_AREA),lumen_centerY(max_AREA),'b.')
% plot(thrombus_centerX(max_AREA),thrombus_centerY(max_AREA),'r.')
% hold off

%Finds 55% mark along the vesse to find centroid data
cent_slice = z_length*.55;
cent_slice = round(cent_slice);
MAG = centroid_diff(cent_slice,1);
ANG = centroid_diff(cent_slice,2);
LEN = cent_slice*voxel_size;

datatable4 = [{'Cent Diff'} MAG;{'Angle'} ANG;{'Length'} LEN]

imshow(thrombus_mesh(:,:,cent_slice))
hold on
plot(lumen_centerX(cent_slice),lumen_centerY(cent_slice),'b.')
plot(thrombus_centerX(cent_slice),thrombus_centerY(cent_slice),'r.')
hold off

%Flips the matrix to start the plot at the proximal point
centroid_diff = flip(centroid_diff);

%plots the magnitude and angle between centroids
% plot(Q1_lumen_area(:,3),centroid_diff(:,1),'.r');
% hold on
% plot(Q1_lumen_area(:,3),centroid_diff(:,2),'.b');
% xlabel('Vessel Length (mm)')
% title('Magnitude and Angle Between Lumen and Outer Lumen Centroids')
% hold off
% legend('Centroid Distance (mm)', 'Angle [0,2pi]')   

% 
% [max_Centdiff,Cent_Dist] = max(centroid_diff(:,1));
% max_Centdiff = Cent_Dist*voxel_size;

max_centdiff = centroid_diff(max_thrombus,1);
max_ang = centroid_diff(max_thrombus,2);

% datatable5 = [{'Cent Diff'} max_centdiff;{'Angle'} max_ang]



%Calculates the percent thrombus in the lumen at each z-slice
% area_percent = zeros(z_length,5);
% for f = 1:z_length   
%    if vessel_area(f) ~=0
%        area_percent(f,1) = (total_thrombusA(f)/vessel_area(f))*100;
%        area_percent(f,2) = (Q1_thrombusA(f)/vessel_area(f))*100;
%        area_percent(f,3) = (Q2_thrombusA(f)/vessel_area(f))*100;
%        area_percent(f,4) = (Q3_thrombusA(f)/vessel_area(f))*100;
%        area_percent(f,5) = (Q4_thrombusA(f)/vessel_area(f))*100;
%    end
% end

% plot(Q1_lumen_area(:,3),area_percent(:,1),'.r');
% hold on
% plot(Q1_lumen_area(:,3),area_percent(:,2),'.b')
% plot(Q1_lumen_area(:,3),area_percent(:,3),'.g')
% plot(Q1_lumen_area(:,3),area_percent(:,4),'.k')
% plot(Q1_lumen_area(:,3),area_percent(:,5),'.m')
% hold off




    
    
    