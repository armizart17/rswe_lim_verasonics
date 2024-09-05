function [grad_z,grad_x,k,sws_matrix] = phase_estimator_lsq(u, w_kernel,f_v,dinf,og_size,constant)
% function [grad_z,grad_x,k,sws_matrix] = phase_estimator_lsq(u, w_kernel,f_v,dinf,og_size,constant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the shear wave speed of a region with the 
% phase gradient method with LEAST SQUARES. 
% 
% Inputs:  
%          u           : 2D region of interest to evaluate (previously mirror padding)
%          w_kernel    : vector with the size of the window kernel
%          f_v         : vibration frequency
%          dinf        : structure that contains the spatial resolutions
%          og_size     : vector containing the original size of the data
%          constant    : constant from equations (0.33 gives good results)

% Outputs: 
%          grad_z       : Gradient matrix for the axial direction
%          grad_x       : Gradient matrix for the lateral direction
%          k            : Total Wavenumber matrix 
%          sws_matrix   : Shear wave speed matrix 
% Author: EAMZ, based on SMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resolution
    
     % Axis of the kernels
    z_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(1))*res_z; 
    x_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(2))*res_x; 
    
    %% Initializing vectors and matrixes
    grad_z = zeros(og_size); 
    grad_x = grad_z; % matrixes containing the estimated wavenumber along each direction
    
    angle_u = angle(u);
    %angle_u_unwrp = unwrap(angle_u, [], direction);
    %angle_u_unwrp = unwrap(angle_u, [], 2);   
    

    [X, Z] = meshgrid(x_axis,z_axis);
    A = [X(:) Z(:) ones(length(x_axis)*length(z_axis),1)];
    At = A';
    AtA = At*A;
    pseudoInverseA = inv(AtA)*At; %LSM

    angle_z = unwrap(angle_u,[],1);
    angle_x = unwrap(angle_u,[],2);
    for ii = 1:og_size(1)
        for jj=1:og_size(2) %% for faster computing pararell toolbox
            area_z = angle_z(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            area_x = angle_x(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            %% Least squares
            b_z = area_z(:);
            b_x = area_x(:);
            results_z = pseudoInverseA * b_z;
            results_x = pseudoInverseA * b_x;
            grad_z(ii,jj) = results_z(2);
            grad_x(ii,jj) = results_x(1);
        end
    end
    
    phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;

    med_wind = floor (2.5/f_v/dinf.dx)*2+1; %the median window contains at least a wavelenght
    k2_med = medfilt2(phase_grad_2,[med_wind med_wind],'symmetric');
    k = sqrt(k2_med);

    sws_matrix = (2*pi*f_v)./k;
end