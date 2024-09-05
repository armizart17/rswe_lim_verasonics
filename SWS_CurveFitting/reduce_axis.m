function [z_axis2,z_vector2,x_axis2,x_vector2] = reduce_axis(track,z_axis,z_vector,x_axis,x_vector)

    location1 = ( abs(z_axis+track) == min(abs(z_axis+track)) );  
    loc1(1) = find(location1==1);
    location1 = ( abs(z_axis-track) == min(abs(z_axis-track)) );  
    loc1(2) = find(location1==1);

    z_axis2 = z_axis(loc1(1):loc1(2));
    z_vector2 = z_vector(loc1(1):loc1(2));


    location1 = (abs(x_axis+track)==min(abs(x_axis+track)));  
    loc1(1) = find(location1==1);
    location1 = (abs(x_axis-track)==min(abs(x_axis-track)));  
    loc1(2) = find(location1==1);

    x_axis2 = x_axis(loc1(1):loc1(2));
    x_vector2 = x_vector(loc1(1):loc1(2));

end