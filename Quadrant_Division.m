function [centroid_line,vert_Q1,vert_Q2,vert_Q3,vert_Q4,rad_Q1,rad_Q2,rad_Q3,rad_Q4,area_Q1,area_Q2,area_Q3,area_Q4] = Quadrant_Division(stl_vertices,voxel_size)


%Sorts the matrix in order of descentind z-values
vert_order = sortrows(stl_vertices,3);

%Find minimum and maximum z
z_min = min(vert_order(:,3));
z_max = max(vert_order(:,3));
z_interval = z_min+(3*voxel_size);

%Initialize a zero matrix to fill with centroids
vert_length = length(vert_order);
z_center = zeros(vert_length,2);
c = 1;

%Centroid calculation
while z_min<=z_max
    x_vals = [];
    y_vals = [];
    a = find(vert_order(:,3)>=z_min & vert_order(:,3)<z_interval); %Finds the elements within one voxel_size
    if isempty(a) == 0 %Determines if the vector is empty before running the for loop
        for b = min(a):max(a) %Creates new vectors with the values from the thrombus_order matrix that are in the range
        x_vals(b) = vert_order(b,1);
        y_vals(b) = vert_order(b,2);
        end
        x_vals = x_vals(x_vals ~=0); %Removes zeros from the vector
        y_vals = y_vals(y_vals ~=0);
        y_midpoint = range(y_vals)/2+min(y_vals);
        x_midpoint = range(x_vals)/2+min(x_vals);
        %Creates a shape that can be readable by the centroid function
        xy_vals = [min(x_vals) x_midpoint max(x_vals) x_midpoint min(x_vals);y_midpoint min(y_vals) y_midpoint max(y_vals) y_midpoint];
        polyin = polyshape(xy_vals(1,:),xy_vals(2,:));
        [x,y] = centroid(polyin);
        %Fills values into the z_center matrix
        z_center(c,1) = x;
        z_center(c,2) = y;
        c = c+1; 
        z_min = z_min+(3*voxel_size);
        z_interval = z_interval+(3*voxel_size);
    else
        z_min = z_min+(3*voxel_size);
        z_interval = z_interval+(3*voxel_size);
    end
end

%Remove zeros from the z_center matrix
deleteZrows = z_center(:,1)<=0;
z_center(deleteZrows,:) = [];

%Divide the range of z-values based on the length of z_center
plot_range = range(vert_order(:,3));
Int = plot_range/length(z_center);
z_min = min(vert_order(:,3));
z_spacing = zeros(length(z_center),1);

%Add z-component to z_center
for t = 1:length(z_center)
    z_spacing(t,1) = z_min;
    z_min = z_min+Int;
end
centroid_line = horzcat(z_center, z_spacing);

% Create a new variable for the minimum and maximum z-value
Z_min = min(centroid_line(:,3));
Z_max = max(centroid_line(:,3));
z_length = length(centroid_line);

% Initialize the STL model to the range of the centerline
vertice_range = find(stl_vertices(:,3)>=Z_min & stl_vertices(:,3)<=Z_max);
vert_z = stl_vertices(vertice_range,:);
vert_z = sortrows(vert_z,3);

% Create empty matrices for each quadrant
vert_Q1 = zeros(length(vertice_range),3);
vert_Q2 = zeros(length(vertice_range),3);
vert_Q3 = zeros(length(vertice_range),3);
vert_Q4 = zeros(length(vertice_range),3);

%fill the lumen quadrant matrices
for k = 1:(z_length-1)
    %finds the elements in lumen_z that are in the range of z-values
    index = find(vert_z(:,3)>=centroid_line(k,3) & vert_z(:,3)<centroid_line(k+1,3));
    for d = min(index):max(index) %Sorts through the elements of lumen_z to pick out where the lie around the centerline
        if vert_z(d,1)<centroid_line(k,1) && vert_z(d,2)<centroid_line(k,2)
            vert_Q2(d,:) = vert_z(d,:);
        elseif vert_z(d,1)<centroid_line(k,1) && vert_z(d,2)>=centroid_line(k,2)
            vert_Q4(d,:) = vert_z(d,:);
        elseif vert_z(d,1)>=centroid_line(k,1) && vert_z(d,2)<centroid_line(k,2)
            vert_Q1(d,:) = vert_z(d,:);
        else vert_z(d,1)>=centroid_line(k,1) && vert_z(d,2)>=centroid_line(k,2);
            vert_Q3(d,:) = vert_z(d,:);
        end
    end
end

% Deletes the zeros out of each matrix
deleteQ1rows = vert_Q1(:,1)==0;
vert_Q1(deleteQ1rows,:) = [];
deleteQ2rows = vert_Q2(:,1)==0;
vert_Q2(deleteQ2rows,:) = [];
deleteQ3rows = vert_Q3(:,1)==0;
vert_Q3(deleteQ3rows,:) = [];
deleteQ4rows = vert_Q4(:,1)==0;
vert_Q4(deleteQ4rows,:) = [];

%Quadrant radius vectors
rad_Q1 = [];
rad_Q2 = [];
rad_Q3 = [];
rad_Q4 = [];

%calculates the radius of each z-slice along the quadrants
for s = 1:(z_length-1)
    %Initialize the positions to calculate lumen radii
    Q1_dist = [];
    Q1_int = find(vert_Q1(:,3)>=centroid_line(s,3) & vert_Q1(:,3)<centroid_line(s+1,3));
    Q2_dist = [];
    Q2_int = find(vert_Q2(:,3)>=centroid_line(s,3) & vert_Q2(:,3)<centroid_line(s+1,3));
    Q3_dist = [];
    Q3_int = find(vert_Q3(:,3)>=centroid_line(s,3) & vert_Q3(:,3)<centroid_line(s+1,3));
    Q4_dist = [];
    Q4_int = find(vert_Q4(:,3)>=centroid_line(s,3) & vert_Q4(:,3)<centroid_line(s+1,3));

    %Finds radii for Q1
    for j = min(Q1_int):max(Q1_int)
        x_Q1 = abs(vert_Q1(j,1))-centroid_line(s,1);
        y_Q1 = abs(vert_Q1(j,2))-centroid_line(s,2);
        rad_Q1 = sqrt(((x_Q1^2)+(y_Q1^2)));
        Q1_dist(j) = rad_Q1;
    end
    Q1_dist = Q1_dist(Q1_dist ~=0);
    rad_Q1(s) = mean(Q1_dist);
    %Finds radii for Q2
    for g = min(Q2_int):max(Q2_int)
        x_Q2 = abs(vert_Q2(g,1))-centroid_line(s,1);
        y_Q2 = abs(vert_Q2(g,2))-centroid_line(s,2);
        rad_Q2 = sqrt(((x_Q2^2)+(y_Q2^2)));
        Q2_dist(g) = rad_Q2;
    end
    Q2_dist = Q2_dist(Q2_dist ~=0);
    rad_Q2(s) = mean(Q2_dist);
    %Finds radii for Q3
    for h = min(Q3_int):max(Q3_int)
        x_Q3 = abs(vert_Q3(h,1))-centroid_line(s,1);
        y_Q3 = abs(vert_Q3(h,2))-centroid_line(s,2);
        rad_Q3 = sqrt(((x_Q3^2)+(y_Q3^2)));
        Q3_dist(h) = rad_Q3;
    end
    Q3_dist = Q3_dist(Q3_dist ~=0);
    rad_Q3(s) = mean(Q3_dist);
    %Finds radii for lumen_Q4
    for r = min(Q4_int):max(Q4_int)
        x_Q4 = abs(vert_Q4(r,1))-centroid_line(s,1);
        y_Q4 = abs(vert_Q4(r,2))-centroid_line(s,2);
        rad_Q4 = sqrt(((x_Q4^2)+(y_Q4^2)));
        Q4_dist(r) = rad_Q4;
    end
    Q4_dist = Q4_dist(Q4_dist ~=0);
    rad_Q4(s) = mean(Q4_dist);

end

%Calculate quadrant area assuming each is a quater circle
area_Q1 = (pi/4)*(rad_Q1.^2);
area_Q2 = (pi/4)*(rad_Q2.^2);
area_Q3 = (pi/4)*(rad_Q3.^2);
area_Q4 = (pi/4)*(rad_Q4.^2);


end




