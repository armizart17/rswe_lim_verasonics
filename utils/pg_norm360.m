function [grad_abs, size_out, grad_avg] = pg_norm360(u, w_kernel, dinf, og_size, stride, num_angles)
% function [grad_abs, size_out] = pg_norm360(u, w_kernel, dinf, og_size, stride)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the L2-norm of the phase gradient for a denoising
% approach considering a 360º rotation for angular average in Kx and Kz
% From the original phase gradient framework by J. Ormachea et al.
% Considering: || \grad(\phi)) ||_2 = || k||_2, -----> eq. \Phi = K
% Later optimize (wavenumber K): min || K - \Phi ||^2_2 + R(K)
% with R (regularization prior, i.e. TV)
% 
% Inputs:  
%          u           : 2D region of interest to evaluate (previously mirror padding)
%          w_kernel    : vector with the size of the window kernel
%          f_v         : vibration frequency
%          dinf        : structure that contains the spatial resolutions
%          og_size     : vector containing the original size of the data
%          stride      : stride for window
%          num_angles  : number of angles for 360º rotation
% Outputs: 
%          grad_abs    : 
%          size_out    : size of grad_abs considering stride and kernel size
% Author: EMZ + GFB
% Explanation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % Function to calculate L2-norm of phase gradient using a given number of angles,
    % handling potential decimal issues with delta angles.

    % Inputs parsing and initial setup
    st = stride;
    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resolution

    % Axis of the kernels
    z_axis = linspace(-(w_kernel(1)-1)/2, (w_kernel(1)-1)/2, w_kernel(1)) * res_z;
    x_axis = linspace(-(w_kernel(2)-1)/2, (w_kernel(2)-1)/2, w_kernel(2)) * res_x;
    
    % Initializing phase angles and output size
    angle_u = angle(u);
    size_mirror = size(u);
    numkernels = floor((size_mirror - w_kernel) ./ st + 1);
    size_out = floor((og_size - 1) ./ st + 1);

    angle_z = unwrap(angle_u, [], 1);
    angle_x = unwrap(angle_u, [], 2);
    
    angle_2d = unwrap_phase(angle_u); % v2

    % figure, 
    % subplot(2,2,1), imagesc(angle_u); title('angle u')
    % subplot(2,2,3), imagesc(angle_z); title('angle z')
    % subplot(2,2,4), imagesc(angle_x); title('angle x')
    % subplot(2,2,2), imagesc(angle_2d); title('angle 2d')

    grad_abs = zeros(size_out);

    % Compute the search space once
    search_z = 1:st:og_size(1);
    search_x = 1:st:og_size(2);

    % Calculate the rotation angle step based on the number of angles
    delta_theta = 360 / num_angles;  % This might be a decimal value

    % Pre-allocate 3D array for storing gradients across rotations
    val_abs_3D = zeros(w_kernel(1), w_kernel(2), num_angles);

    % Main loop for gradient calculation
    for ii = 1:length(search_z)
        for jj = 1:length(search_x)
            
            % Extract v1
            % area_z = angle_z(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1); 
            % area_x = angle_x(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1);

            % Extract v2
            % area_z = angle_2d(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1); 
            % area_x = area_z;
             
            % if(ii == fix( length(search_z) / 2) && jj == fix( length(search_x) / 2) )
            %     keyboard;
            % end

            % Extract v3 (unwrapp per window)       
            % area_z = angle_u(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1); 
            % area_x = angle_u(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1);
            % 
            % area_z = unwrap(area_z, [], 1);
            % area_x = unwrap(area_x, [], 2);

            % Extract v4 (unwrapp per window 2d) 
            area_2d = angle_u(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1); 

            area_z = unwrap_phase(area_2d);
            area_x = area_z;

            % Parallel loop over rotations using the given number of angles
            for idx = 1:num_angles

                clear diff_area_z diff_area_x
                % Calculate the actual rotation angle
                theta = mod((idx - 1) * delta_theta, 360); % Ensuring the angle wraps around at 360

                % Rotate areas and calculate gradients
                rotated_area_z = imrotate(area_z, theta, 'nearest', 'crop');
                rotated_area_x = imrotate(area_x, theta, 'nearest', 'crop');
                % 
                % figure, 
                % subplot(221), imagesc(area_x), title('x'), colorbar
                % subplot(222), imagesc(area_z), title('z'), colorbar
                % subplot(223), imagesc(rotated_area_x), title('Rot x'), colorbar
                % subplot(224), imagesc(rotated_area_z), title('Rot z'), colorbar

                diff_area_z = diff(rotated_area_z, 1, 1) ./ res_z;
                diff_area_x = diff(rotated_area_x, 1, 2) ./ res_x;

                % Ensure consistent dimensions by padding the last row/column
                diff_area_z = [diff_area_z; diff_area_z(end, :)]; % replicate last row
                diff_area_x = [diff_area_x, diff_area_x(:, end)];  % replicate last column

                % Compute the L2-norm of the gradients
                % val_abs_3D(:, :, idx) = sqrt(diff_area_z.^2 + diff_area_x.^2);
                aux = sqrt(diff_area_z.^2 + diff_area_x.^2);

                val_abs_3D(:, :, idx) = medfilt2(aux,[7 7],'symmetric');
                % val_abs_3D(:, :, idx) = aux;   % *  
                % val_avg_3D(:, :, idx) = 0.5*(diff_area_z + diff_area_x);
            end

            % Average L2-norm across all rotations
            val = mean(val_abs_3D, 3);
            % val2 = mean(val_avg_3D, 3);

            % grad_abs(ii, jj) = mean(val, 'all'); % *
            grad_abs(ii, jj) = median(val, 'all');
            % grad_avg(ii, jj) = mean(val2, 'all');
        end
    end
end
