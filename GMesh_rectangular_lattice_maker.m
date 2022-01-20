%% Initialising the program
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%INPUT USER VALUES HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimension_a = 0.6;    % (the distance from the centre of one circle to the next)
dimension_r = 0.174;  % (the distance from the centre of a circle to the edge of that circle)
grid_width = 4;       % (number of circles in the X direction)
grid_height = 20;      % (number of circles in the Y direction)

output_filename = "Github Test";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Developed by Amy Thomas%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Making the circle plotting variables
circle_num = grid_height*(grid_width+2);
circle_data = zeros(circle_num, 3);

if grid_height == grid_width
    min_dimension = grid_height;
elseif grid_height < grid_width
    min_dimension = grid_height;
else
    min_dimension = grid_width;
end

lattice_height = (grid_height+1) * dimension_a;
lattice_width = (grid_width+1) * dimension_a;

half_dimension_a = dimension_a / 2;
hori_square_num = 2*grid_width+2;
vert_square_num = 2*grid_height+2;

%Making the geometry variables
counter = 0;
point_count = 0;
ref_point_count = 0;
radial_refs = [-1 1; -1 0; -1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1];
corner_matrix = [0 lattice_height 2 1 1; lattice_width lattice_height 1 1 1; lattice_width 0 2 -1 2; 0 0 1 -1 2; 0 0 2 0 1];
corner_number = [];

line_number = 0;
ext_curve_matrix = [];
int_circle_data = [];
int_curve_matrix = [];
loop_check = 0;
side_list = ["Right" "Bottom" "Left" "Top"];

%Making the surface plane variables
square_plane_matrix = zeros(2*vert_square_num+1, 2*hori_square_num+1);
square_count = 0;
halberd_plane_matrix = zeros(grid_height*3, 4+grid_width*3);
halberd_counter = 1;
circle_plane_matrix = zeros(4*grid_height, 2+2*grid_width);

curve_loop_matrix = [];
loop_number = 0;
loop_data = '';
physical_curve_matrix = [];
char_points = [];

%Making the data for each circle
for n = 1:1:grid_height
    for m = 1:1:grid_width+1
        counter = counter + 1;
        circle_data(counter, 1) = (m-1) * (dimension_a);
        circle_data(counter, 2) = (n) * (dimension_a);
        circle_data(counter, 3) = dimension_r;
        
        if m ~= 1
            int_circle_data = [int_circle_data; circle_data(counter, 1) circle_data(counter, 2) circle_data(counter, 3)];
        end  
    end
    counter = counter + 1;
    circle_data(counter, 1) = lattice_width;
    circle_data(counter, 2) = lattice_height - n * (dimension_a);
    circle_data(counter, 3) = dimension_r;
end

%Create .geo file
filename=append(output_filename, '.geo');
output = fopen(filename, 'wt');
fprintf(output, '//GMesh .geo file created using MATLAB code written by Amy Thomas\n');
fprintf(output, 'SetFactory("OpenCASCADE"); \n \n');

%making the reference matrix
reference_points = zeros(vert_square_num+1,  hori_square_num+1);
corner_refs = [ 1 1; 1 hori_square_num+1; vert_square_num+1 hori_square_num+1; vert_square_num+1 1; 1 1];

%working out the distance to split the lines
separation = half_dimension_a - dimension_r;

%working out the diagonal distance from the intersections to the circle
diag_dist = sqrt(2) * half_dimension_a - dimension_r;
lin_dist = diag_dist / sqrt(2);
diag_dist2 = half_dimension_a - lin_dist;

%% Defining the geometry of the points

%Defining the points for the external edges
fprintf(output, '//Defining the points \n \n');

%Starting with the top corner of the right edge and moving in a clockwise
%direction, this loop calculates all the points that make up the external
%edge of the lattice and assigns them to the relevant matrices for later
for m = 1:1:4
    q=0;

    %adding the corner and recording what its number is for later
    point_count = point_count + 1;
    fprintf(output, 'Point(%i) = {%f, %f, 0, 1.0};\n', point_count, corner_matrix(m, 1), corner_matrix(m, 2));
    
    reference_points(corner_refs(m, 1), corner_refs(m, 2)) = point_count;    
    ext_curve_matrix = [ext_curve_matrix; point_count point_count+1 0 m];
    corner_number = [corner_number; point_count];
    char_points = [char_points; point_count];

    %doing the loop for the right and left sides
    if m == 2 || m == 4
 
        %for each circle in the spreadsheet
        for n = 1:1:circle_num

            %checking if the radius crosses outside the lattice for a given
        	%side. This is just a calculation of the absolute value of the
        	%radius+centre coordinates and checking if they exceed 1.      
            if circle_data(n, corner_matrix(m, 3)) + corner_matrix(m, 4) * dimension_r >= lattice_width || circle_data(n, corner_matrix(m, 3)) + corner_matrix(m, 4) * dimension_r <= 0
                q = q + 1;
                
                %finding the intersection coordinates
                coord1 = circle_data(n, 2) + corner_matrix(m-1, 4)*dimension_r;
                coord2 = circle_data(n, 2) - corner_matrix(m-1, 4)*dimension_r;
                coord_mat = [coord1 coord2];

                %adding 7 points
                point_count = point_count + 4;

                %adding the separator point
                fprintf(output, 'Point(%i) = {%f, %f, 0, 1.0};\n', point_count-3, circle_data(n, 1), circle_data(n, 2)+corner_matrix(m-1, 4)*half_dimension_a);
                %adding the first corner point
                fprintf(output, 'Point(%i) = {%f, %f, 0, 1.0};\n', point_count-2, circle_data(n, 1), coord1);
                %adding the first diagonal point
                fprintf(output, 'Point(%i) = {%0.8f, %0.8f, 0, 1.0};\n', point_count-1*corner_matrix(m-1, 4), circle_data(n, 1)-corner_matrix(m-1, 4)*diag_dist2, circle_data(n, 2)+diag_dist2);
                %adding the radial point
                fprintf(output, 'Point(%i) = {%f, %f, 0, 1.0};\n', point_count, circle_data(n, 1)-corner_matrix(m-1, 4)*dimension_r, circle_data(n, 2));
                %adding the second diagonal point
                fprintf(output, 'Point(%i) = {%0.8f, %0.8f, 0, 1.0};\n', point_count+1*corner_matrix(m-1, 4), circle_data(n, 1)-corner_matrix(m-1, 4)*diag_dist2, circle_data(n, 2)-diag_dist2);            
                %adding the second corner point
                fprintf(output, 'Point(%i) = {%f, %f, 0, 1.0};\n', point_count+2, corner_matrix(m, 1), coord2);
                %adding the centre point            
                fprintf(output, 'Point(%i) = {%f, %f, 0, 1.0};\n', point_count+3, circle_data(n, 1), circle_data(n, 2));            

                %adding all the points to the external curve matrix
                ext_curve_matrix = [ext_curve_matrix; point_count-3 point_count-2 0 m];
                
                %adding the relevant points to the reference matrix
                reference_points(corner_refs(m,1)+corner_matrix(m-1, 4)*(2*q-1), corner_refs(m,2)) = point_count-3;
                reference_points(corner_refs(m,1)+corner_matrix(m-1, 4)*(2*q), corner_refs(m,2)) = point_count+3;
                
                for o = 1:1:4
                    ext_curve_matrix = [ext_curve_matrix; point_count-3+o point_count-2+o point_count+3 m];
                end

                ext_curve_matrix = [ext_curve_matrix; point_count+2 point_count+4 0 m];

                %adding the two corner points to the characteristic points
                %matrix
                char_points = [char_points; point_count-2*corner_matrix(m, 1); point_count+2*corner_matrix(m, 1)];
                
                point_count = point_count + 3;
                
            end
        end
        
        point_count = point_count + 1;
        fprintf(output, 'Point(%i) = {%f, %f, 0, 1.0};\n', point_count, corner_matrix(m, 1), corner_matrix(m, 2)-corner_matrix(m, 4)*half_dimension_a*(vert_square_num-1));
        
        if m == 4 && corner_refs(m+1,1)-corner_matrix(m, 4) == 2
            ext_curve_matrix = [ext_curve_matrix; point_count 1 0 m];
            reference_points(corner_refs(m+1,1)-corner_matrix(m-1, 4), corner_refs(m,2)) = point_count;
        else
            ext_curve_matrix = [ext_curve_matrix; point_count point_count+1 0 m];
            reference_points(corner_refs(m+1,1)-corner_matrix(m-1, 4), corner_refs(m,2)) = point_count;
        end
    end
    
    %doing the loop for the bottom and top sides
    if m == 1 || m ==3                
        for n = 1:1:hori_square_num-1
                    
            point_count = point_count + 1;
            fprintf(output, 'Point(%i) = {%f, %f, 0, 1.0};\n', point_count, corner_matrix(m, 1) + corner_matrix(m, 4) * n * half_dimension_a, corner_matrix(m, 2));
            
            reference_points(corner_refs(m, 1), corner_refs(m, 2)+corner_matrix(m, 4)*n) = point_count;
            ext_curve_matrix = [ext_curve_matrix; point_count point_count+1 0 m];
        end        
    end       
end

%Starting in the top right corner moving downwards this loop calculates
%where all the reference points need to go and links them up to the correct
%other points for later

for n = 1:1:hori_square_num/2
    
    x_point = (2*n-1) * half_dimension_a;
    
    for m = 1:1:vert_square_num-1
        
           y_point = m*half_dimension_a;
           
           point_count = point_count + 1;
           fprintf(output, 'Point(%i) = {%f, %f, 0, 1.0};\n', point_count, x_point, y_point);
           reference_points(vert_square_num+1-m , 2*n) = point_count;

    end
end
for n = 1:1:grid_width
    
    x_point = n * 2 * half_dimension_a;
    
    for m = 1:1:grid_height+1
        
           y_point = (2*m-1)*half_dimension_a;
           
           point_count = point_count + 1;
           fprintf(output, 'Point(%i) = {%f, %f, 0, 1.0};\n', point_count, x_point, y_point);
           reference_points(vert_square_num+2-2*m , 1+2*n) = point_count;

    end    
end

%Starting with the circle in the top right and moving in a clockwise
%direction, this loop calculates all the points that make up the internal
%circles and assigns them to the relevant matrices for later

int_size = size(int_circle_data);
for n = 1:1:int_size(1,1)
    
    %printing the centre point of each internal circle
    point_count = point_count + 1;
    fprintf(output, 'Point(%i) = {%0.5f, %0.5f, 0, 1.0};\n', point_count, int_circle_data(n, 1), int_circle_data(n, 2));
    
    z = 0;
    for o = vert_square_num+1:-1:1
    	for q = 1:1:hori_square_num+1
            if reference_points(o, q) == 0
                reference_points(o, q) = point_count;
                z = 1;
                break
            end
        end
        if z == 1
        	break
        end
    end
    
    for m = 1:1:8
        
        fprintf(output, 'Point(%i) = {%0.8f, %0.8f, 0, 1.0};\n', point_count+m, int_circle_data(n, 1)+dimension_r*cos(m*pi/4),  int_circle_data(n, 2)+dimension_r*sin(m*pi/4));
        if m == 8
            int_curve_matrix = [int_curve_matrix; point_count+m point_count+1 point_count];
        else
            int_curve_matrix = [int_curve_matrix; point_count+m point_count+m+1 point_count];            
        end
    end
    
    point_count = point_count + 8; 
    
end
corner_number = [corner_number; point_count];

%% Defining the geometry of the lines
% Connecting the appropriate points to draw the lattice shape and reference
% lines correctly

fprintf(output, '\n//Defining the geometric lines and curves \n \n');

%defining the external lines and curves
ext_size = size(ext_curve_matrix);
position = [hori_square_num 1];
circle_count = 0;
halberd_count = 0;
for n = 1:1:ext_size(1,1)
    if ext_curve_matrix(n, 3) == 0

        line_number = line_number + 1;
        fprintf(output, 'Line(%i) = {%i, %i};\n', line_number, ext_curve_matrix(n, 1), ext_curve_matrix(n, 2));
        
        %putting the line into the sqaure matrix
        if ext_curve_matrix(n, 4) == 1
            square_plane_matrix(1, 2*line_number) = line_number;
        elseif ext_curve_matrix(n, 4) == 3
            square_count = square_count + 1;
            square_plane_matrix(1+vert_square_num*2, 2+2*(hori_square_num-square_count)) = line_number;
        elseif ext_curve_matrix(n, 4) == 2
            halberd_count = halberd_count + 1;
            if halberd_count == 1 || halberd_count == 3*grid_height+2
            else
              halberd_plane_matrix(halberd_count-1, 4+grid_width*3) = line_number;
            end
        else
            halberd_count = halberd_count + 1;
            if halberd_count == 3*grid_height+3 || halberd_count == 6*grid_height+4
            else
              halberd_plane_matrix(6*grid_height+4-halberd_count, 1) = line_number;
              
            end
        end
    else
        
        line_number = line_number + 1;
        halberd_count = halberd_count + 0.25;
        fprintf(output, 'Circle(%i) = {%i, %i, %i};\n', line_number,  ext_curve_matrix(n, 1), ext_curve_matrix(n, 3), ext_curve_matrix(n, 2));
        
        if ext_curve_matrix(n, 4)==2
            circle_count = circle_count + 1;
            circle_plane_matrix(circle_count, hori_square_num) = line_number;
        else
            circle_plane_matrix(circle_count, 1) = line_number;
            circle_count = circle_count - 1;            
        end
    end
    
    square_plane_matrix(2, 1+2*hori_square_num) = hori_square_num+1;
    square_plane_matrix(2*vert_square_num, 1+2*hori_square_num) = hori_square_num+2+6*grid_height;   
    square_plane_matrix(2*vert_square_num, 1) = 2*hori_square_num+3+6*grid_height;   
    square_plane_matrix(2, 1) = ext_size(1,1);
end

%defining the internal circle curves
int_size = size(int_curve_matrix);
loop_check1 = line_number;
circle_coord = [0 0;0 -1; 1 -1; 2 -1; 3 -1; 3 0; 2 0 ; 1 0];
circle_count = 0;
hori_count = 1;
vert_count = 1;
for n = 1:1:int_size(1,1)
	
    line_number = line_number + 1;
    circle_count = circle_count + 1;
    fprintf(output, 'Circle(%i) = {%i, %i, %i};\n', line_number,  int_curve_matrix(n, 1), int_curve_matrix(n, 3), int_curve_matrix(n, 2));    
    
    circle_plane_matrix(1+4*grid_height-4*vert_count+circle_coord(circle_count, 1), 1+2*hori_count+circle_coord(circle_count, 2)) = line_number;
    
    if circle_count == 8
        circle_count = 0;
        if hori_count == grid_width
            hori_count = 1;
            vert_count = vert_count + 1;
        else
            hori_count = hori_count + 1;
        end
    end
end

%defining the grid reference lines
for n = 1:1:hori_square_num/2
	for m = 1:1:vert_square_num
    	line_number = line_number + 1;
    	fprintf(output, 'Line(%i) = {%i, %i};\n', line_number, reference_points(m, 2*n), reference_points(m+1, 2*n));
    
        square_plane_matrix(2*m, 2*(2*n)-1) = line_number;
    end   
end 
for n = 1:1:grid_width
    line_number = line_number + 2;
	fprintf(output, 'Line(%i) = {%i, %i};\n', line_number-1, reference_points(1, 2*n+1), reference_points(2, 2*n+1));
	fprintf(output, 'Line(%i) = {%i, %i};\n', line_number, reference_points(vert_square_num, 2*n+1), reference_points(vert_square_num+1, 2*n+1));
   
    square_plane_matrix(2, 4*n+1) = line_number-1;
    square_plane_matrix(2*vert_square_num, 4*n+1) = line_number;
end
for n = 1:1:grid_height+1
    for m = 1:1:hori_square_num
        
        line_number = line_number + 1;
        fprintf(output, 'Line(%i) = {%i, %i};\n', line_number, reference_points(2*n, m), reference_points(2*n, m+1));
        
        square_plane_matrix(3+4*(n-1), 2*m) = line_number;
    end
end 

%defining the radial reference lines
for n = 1:1:grid_width+2
    for m = 1:1:grid_height
        coord1 = 2*m+1;
        coord2 = 2*n-1;       
        centre_point1 = reference_points(coord1, coord2);
        if n == 1
            for o = 1:1:3
                line_number = line_number + 1;
                fprintf(output, 'Line(%i) = {%i, %i};\n', line_number, centre_point1-5+o, reference_points(coord1 + radial_refs(6+o, 1), coord2 + radial_refs(6+o, 2)));
                
                halberd_plane_matrix(3*m-o+1, 2) = line_number;
            end
        elseif n == grid_width+2
            for o = 1:1:3
                line_number = line_number + 1;
                fprintf(output, 'Line(%i) = {%i, %i};\n', line_number, centre_point1-5+o, reference_points(coord1 + radial_refs(2+o, 1), coord2 + radial_refs(2+o, 2)));

                halberd_plane_matrix(3*(m-1)+o, 3+3*grid_width) = line_number;
            end
        else    
            for o = 1:1:8
                line_number = line_number + 1;
                fprintf(output, 'Line(%i) = {%i, %i};\n', line_number, centre_point1+o, reference_points(coord1 + radial_refs(o, 1), coord2 + radial_refs(o, 2)));
                
                halberd_plane_matrix(3*m-1+radial_refs(o, 1), 1+3*(n-1)+radial_refs(o, 2)) = line_number;
            end
        end
    end
end

%% Defining the plane surfaces, curve loops and characteristic geometries
% Drawing the correct data together to define the lattice to allow it to be
% meshed correctly

%Define plane surface and boundary
fprintf(output, '\n//Defining the curve loops and plane surfaces \n \n');

%writing the physical curves for the square and halberd matrices
square_count = [2; 2*vert_square_num];
for m = 1:1:2
    for n = 1:1:hori_square_num
        loop_number = loop_number + 1;
        loop_data = '';
        loop_line = [square_plane_matrix(square_count(m)-1, 2*n) square_plane_matrix(square_count(m), 2*n+1) square_plane_matrix(square_count(m)+1, 2*n) square_plane_matrix(square_count(m), 2*n-1)];

        for o = 1:1:4
            string = strcat(mat2str(loop_line(1, o)),',  ');
            loop_data = strcat(loop_data, string);
        end
        loop_data = loop_data(1:end-1);
        fprintf(output, 'Curve Loop(%i) = {%s};\n', loop_number, loop_data);

    end
end

halberd_coord = [-1 0;-1 1;0 1;1 1; 1 0];
halberd_coord2 = [-1 0;-1 -1;0 -1;1 -1; 1 0];
square_coord = [-2 1;-1 2; 1 2; 2 1];
square_coord2 = [-2 -1;-1 -2; 1 -2; 2 -1];
for n = 1:1:hori_square_num
    for o = 1:1:grid_height
        for q = 1:1:4
            loop_number = loop_number + 1;
            loop_data = '';
            
            if rem(n, 2)==1
                loop_line = [circle_plane_matrix(q+4*(o-1), n) halberd_plane_matrix(2+3*(o-1)+halberd_coord(q, 1), 1+1.5*(n-1)+halberd_coord(q, 2)) square_plane_matrix(5+4*(o-1)+square_coord(q, 1), 1+2*(n-1)+square_coord(q, 2)) halberd_plane_matrix(2+3*(o-1)+halberd_coord(q+1, 1), 1+1.5*(n-1)+halberd_coord(q+1, 2))];
            else
                loop_line = [circle_plane_matrix(q+4*(o-1), n) halberd_plane_matrix(2+3*(o-1)+halberd_coord2(q, 1),4+3*(0.5*n-1)+halberd_coord2(q, 2)) square_plane_matrix(5+4*(o-1)+square_coord2(q, 1), 5+4*(0.5*n-1)+square_coord2(q, 2)) halberd_plane_matrix(2+3*(o-1)+halberd_coord2(q+1, 1), 4+3*(0.5*n-1)+halberd_coord2(q+1, 2))];
            end
            for m = 1:1:4
                string = strcat(mat2str(loop_line(1, m)),',  ');
                loop_data = strcat(loop_data, string);
            end
            loop_data = loop_data(1:end-1);
            fprintf(output, 'Curve Loop(%i) = {%s};\n', loop_number, loop_data);
        end
    end
end

for n = 1:1:loop_number
    fprintf(output, 'Plane Surface(%i) = {%i};\n', n, n);
end

fprintf(output, "For n In {1:%i}\n    Transfinite Surface {n};\n	Recombine Surface {n};\nEndFor", hori_square_num*2);

%finishing up
fclose(output);