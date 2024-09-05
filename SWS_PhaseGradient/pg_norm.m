function [grad_abs, size_out] = pg_norm(u, w_kernel, dinf, og_size, stride)
% function [grad_abs, size_out] = pg_norm(u, w_kernel, dinf, og_size, stride)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the L2-norm of the phase gradient for a denoising approach
% From the original phase gradient framework by J. Ormachea et al.
% Considering: || \grad(\phi)) ||_2 = || k||_2, -----> eq. \Phi = K
% Later optimize (wavenumber K): min || K - \Phi ||^2_2 + R(K)
% with R (regularization prior, i.e. TV)
% PHASE GRADIENT FRAMEWORK BY J. ORMACHEA
% 
% Inputs:  
%          u           : 2D region of interest to evaluate (previously mirror padding)
%          w_kernel    : vector with the size of the window kernel
%          f_v         : vibration frequency
%          dinf        : structure that contains the spatial resolutions
%          og_size     : vector containing the original size of the data
%          stride      : stride for window

% Outputs: 
%          grad_abs    : 
%          size_out    : size of grad_abs considering stride and kernel  size
% Author: EMZ 4 S.Merino
% Explanation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    st = stride;
    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resolution
    
     % Axis of the kernels
    z_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(1))*res_z; 
    x_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(2))*res_x; 
    
    %% Initializing vectors and matrixes    
    angle_u = angle(u);
    
    %% HELP FOR DIMENSIONS %% % THEORY
    size_mirror = size(u); % ogsize + w_kernel - 1; %  = size (u)
    numkernels = floor( (size_mirror - w_kernel)./st + 1 ); % numKernels
    size_out = floor( (og_size - 1)./st + 1 );
   	     
    cont_kernel = 1; 
    angle_z = unwrap(angle_u,[],1);
    angle_x = unwrap(angle_u,[],2);

    grad_abs = zeros(size_out);
    
    %% NORMAL LOOP
    i_out = 1;
    for ii = 1:st:og_size(1)

        j_out = 1;

        for jj = 1:st:og_size(2) %% for faster computing pararell toolbox

            area_z = angle_z(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            area_x = angle_x(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel

            diff_area_z = diff(area_z, 1, 1) ./ res_z;
            diff_area_x = diff(area_x, 1, 2) ./ res_x;

            diff_area_z = [diff_area_z; diff_area_z(end, :)];  % replicate the last row
            diff_area_x = [diff_area_x, diff_area_x(:, end)];  % replicate the last column

            val_abs = sqrt(diff_area_z.^2 + diff_area_x.^2);
            grad_abs(i_out,j_out) = mean(val_abs, 'all');

            j_out = j_out + 1;

        end
       i_out = i_out + 1;
    end

    %% PARALLEL MODIFICATION LOOP
    % search_z = 1 : st : og_size(1);
    % 
    % for ii = 1 : length(search_z)
    % 
    %     search_x = 1 : st : og_size(2);
    % 
    %     parfor jj = 1 : length(search_x)
    % 
    %         area_z = angle_z(search_z(ii):search_z(ii)+w_kernel(1)-1, search_x(jj):search_x(jj)+w_kernel(2)-1); % Window kernel
    %         area_x = angle_x(search_z(ii):search_z(ii)+w_kernel(1)-1, search_x(jj):search_x(jj)+w_kernel(2)-1); % Window kernel
    % 
    %         diff_area_z = diff(area_z, 1, 1) ./ res_z;
    %         diff_area_x = diff(area_x, 1, 2) ./ res_x;
    % 
    %         diff_area_z = [diff_area_z; diff_area_z(end, :)];  % replicate the last row
    %         diff_area_x = [diff_area_x, diff_area_x(:, end)];  % replicate the last column
    % 
    %         val_abs = sqrt(diff_area_z.^2 + diff_area_x.^2);
    %         grad_abs(ii,jj) = mean(val_abs, 'all');    
    % 
    %     end
    % 
    % 
    % end

    
end