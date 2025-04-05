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
    [vert_thrombus_d,faces_thrombus_d,normals_thrombus_d,name_thrombus_d] = stlRead([path thrombus_stl]);

    
%Runs the STL files through the Quadrant_Division function
[lumen_centroid,lumen_Q1,lumen_Q2,lumen_Q3,lumen_Q4,lumen_Q1_rad,lumen_Q2_rad,lumen_Q3_rad,lumen_Q4_rad,lumen_Q1_area,lumen_Q2_area,lumen_Q3_area,lumen_Q4_area] = Quadrant_Division(vert_lumen_d,voxel_size);
[thrombus_centroid,thrombus_Q1,thrombus_Q2,thrombus_Q3,thrombus_Q4,thrombus_Q1_rad,thrombus_Q2_rad,thrombus_Q3_rad,thrombus_Q4_rad,thrombus_Q1_area,thrombus_Q2_area,thrombus_Q3_area,thrombus_Q4_area] = Quadrant_Division(vert_thrombus_d,voxel_size);


% plot3(thrombus_Q1(:,1),thrombus_Q1(:,2),thrombus_Q1(:,3),'r.');
% hold on
% plot3(thrombus_Q2(:,1),thrombus_Q2(:,2),thrombus_Q2(:,3),'g.');
% plot3(thrombus_Q3(:,1),thrombus_Q3(:,2),thrombus_Q3(:,3),'b.');
% plot3(thrombus_Q4(:,1),thrombus_Q4(:,2),thrombus_Q4(:,3),'c.');
% plot3(thrombus_centroid(:,1),thrombus_centroid(:,2),thrombus_centroid(:,3),'k.');
% hold off
    
plot3(lumen_Q1(:,1),lumen_Q1(:,2),lumen_Q1(:,3),'r.');
hold on
plot3(lumen_Q2(:,1),lumen_Q2(:,2),lumen_Q2(:,3),'g.');
plot3(lumen_Q3(:,1),lumen_Q3(:,2),lumen_Q3(:,3),'b.');
plot3(lumen_Q4(:,1),lumen_Q4(:,2),lumen_Q4(:,3),'c.');
plot3(lumen_centroid(:,1),lumen_centroid(:,2),lumen_centroid(:,3),'k.');
hold off   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    