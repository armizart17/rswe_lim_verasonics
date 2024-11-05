function [grad_l2, size_out] = pg_norm_rot(u, w_kernel, dinf, og_size, stride, num_angles, version)
% function [grad_l2, size_out] = pg_norm_rot(u, w_kernel, dinf, og_size, stride, version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the L2-norm of the phase gradient for a denoising
% approach considering a 180º ROTATION for angular average in Kx and Kz defined
% by num_angles. To avoid issues it uses a circle mask of the window
% From the original phase gradient framework by J. Ormachea et al.
% Considering: || \grad(\phi)) ||_2 = || k ||_2, -----> eq. \Phi = K
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
%          num_angles  : number of angles for a 360º rotation
%          version (optional): 
%                   1 Unwrapp XZ ALL angle part. velocity
%                   2 Unwrapp 2D function ALL angle part. velocity
%                   3 Unwrapp XZ each WINDOW angle part. velocity (default)
%                   4 Unwrapp 2D each WINDOW angle part. velocity
% Outputs: 
%          grad_l2     : grad eq K 
%          size_out    : size of grad_l2 considering stride and kernel size
% Author: EMZ 
% Notes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % Function to calculate L2-norm of phase gradient using a given number of angles,
    % handling potential decimal issues with delta angles.

    if nargin < 7  % If less than 7 arguments, set 'version' to 3
        version = 3;
    end


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

    angle_z = unwrap(angle_u, [], 1); % v1
    angle_x = unwrap(angle_u, [], 2); % v1
    
    angle_2d = unwrap_phase(angle_u); % v2

    grad_l2 = zeros(size_out);

    % Compute the search space once
    search_z = 1:st:og_size(1);
    search_x = 1:st:og_size(2);

    % Calculate the rotation angle step based on the number of angles
    delta_theta = 360 / num_angles;  % This might be a decimal value

    % Pre-allocate 3D array for storing gradients across rotations
    val_abs_3D = zeros(w_kernel(1), w_kernel(2), num_angles);

    % Circle mask for rotation 
    M = w_kernel(1); 
    N = w_kernel(2);
    [xGrid, yGrid] = meshgrid(1:N, 1:M);
    
    % Check the number of angles to determine whether to apply a circle mask or full ones
    if ismember(num_angles, [1, 2, 4]) % (90º multiple angles)
        % Use a full mask of ones if num_angles is 1, 2, 
        circ_mask = logical(ones(M, N));
    else
        % Calculate the circle mask for other values of num_angles
        R = floor(min(M, N) / 2);  % Radius of the circle
        centerX = ceil(N / 2);     % X-coordinate of the center
        centerY = ceil(M / 2);     % Y-coordinate of the center
        
        % Create the circular mask
        circ_mask = sqrt((xGrid - centerX).^2 + (yGrid - centerY).^2) < R;
    end

    % Main loop for gradient calculation
    for ii = 1:length(search_z)
        for jj = 1:length(search_x)
            
            switch version
                case 1
                    % Extract v1
                    area_z = angle_z(search_z(ii):search_z(ii)+w_kernel(1)-1, search_x(jj):search_x(jj)+w_kernel(2)-1); 
                    area_x = angle_x(search_z(ii):search_z(ii)+w_kernel(1)-1, search_x(jj):search_x(jj)+w_kernel(2)-1);
            
                case 2
                    % Extract v2
                    area_z = angle_2d(search_z(ii):search_z(ii)+w_kernel(1)-1, search_x(jj):search_x(jj)+w_kernel(2)-1); 
                    area_x = area_z;
            
                case 3
                    % Extract v3 (unwrap per window) 
                    area_z = angle_u(search_z(ii):search_z(ii)+w_kernel(1)-1, search_x(jj):search_x(jj)+w_kernel(2)-1); 
                    area_x = angle_u(search_z(ii):search_z(ii)+w_kernel(1)-1, search_x(jj):search_x(jj)+w_kernel(2)-1);
                    
                    % Unwrapping
                    area_z = unwrap(area_z, [], 1);  % Unwrap in dimension 1
                    area_x = unwrap(area_x, [], 2);  % Unwrap in dimension 2
            
                case 4
                    % Extract v4 (unwrap per window 2D)
                    area_2d = angle_u(search_z(ii):search_z(ii)+w_kernel(1)-1, search_x(jj):search_x(jj)+w_kernel(2)-1); 
                    
                    % Use a custom phase unwrapping function
                    area_z = unwrap_phase(area_2d);
                    area_x = area_z;
            end

            % Parallel loop over rotations using the given number of angles
            for idx = 1:num_angles
            % parfor idx = 1:num_angles % faster**

                clear diff_area_z diff_area_x
                % Calculate the actual rotation angle
                theta = mod((idx - 1) * delta_theta, 180); % Ensuring the angle wraps around at 360

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
                val_abs_3D(:, :, idx) = sqrt(diff_area_z.^2 + diff_area_x.^2);

     
            end

            % Average L2-norm across all rotations through depth (angle dimension)
            val = mean(val_abs_3D, 3); % 2D-array

            val_circ_mask = val(circ_mask); % select circle mask returns a 1D-array

            grad_l2(ii,jj) = mean(val_circ_mask,'all');
        end
    end
end
