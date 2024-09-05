function sws_matrix = sws_generator(u,window,f_v,d_n,dinf,og_size,est_z,est_x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functition that yields the shear wave speed of a region with the 
% approximation method. 
% 
% 
% Inputs:  u           - 2D region of interest to evaluate
%          window      - vector with the size of the window kernel
%          f_v         - vibration freqeuncy %700 
%          d_n         - Number of pixels to define the lag
%          dinf        - structure that contains the spatial resolutions
%          og_size     - vector containing the original size of the data
%          est_z       - Axial estimator (usually 10)
%          est_x       - Lateral estimator (usually 5)
%
% Outputs: sws_matrix  - Resulting shear wave speed matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resoltion
    delta_z = d_n*res_z; % Axial lag
    delta_x = d_n*res_x; % Lateral lag
    c_z = est_z/(delta_z^2); % Axial constant 
    c_x = est_x/(delta_x^2); % Lateral constant 
    T_A = xcorr2(ones(window)); % Correction triangle

    % Initializations of the wavenumber and SWS matrixes
    k_z = zeros(length(1:2:og_size(1)),length(1:2:og_size(2)));
    k_x = k_z;
    sws_matrix = k_x;
    step = 2;

    %% Begin

    search_z = 1:step:og_size(1); % To save time, the windowed matrix is evaluated in points separated by a step of 2
    for ii = 1:length(search_z)

        search_x = 1:step:og_size(2);
        parfor jj=1:length(search_x)
        % for jj=1:length(search_x)
            %area = u(ii:ii+window(1)-1,jj:jj+window(2)-1);
            area = u(search_z(ii):search_z(ii)+window(1)-1,search_x(jj):search_x(jj)+window(2)-1); % Window kernel 
            auto_matrix = xcorr2(area);
            auto_matrix = auto_matrix./T_A; % Discrete function correction
            %center = auto_matrix(window(1),window(2));

            z_vector = auto_matrix(:,window(2))'; % Axial vector from the center of the autocorrelation kernel
            z_vector = z_vector/max(real(z_vector)); % WE CHOOSE MAX AMZ
%             z_vector = z_vector/real(z_vector(window(1)));
            %center_z = z_vector(window(1));
            z_lag = z_vector(window(1) - d_n);
            %k_z(ii,jj) = sqrt( c_z*( 1 - real(z_lag)/real(center_z) ) ); % Axial wavenumber matrix calculation with the approximation
            k_z(ii,jj) = sqrt( c_z*( 1 - real(z_lag) ) ); % Axial wavenumber matrix calculation with the approximation

            x_vector = auto_matrix(window(1),:); % Lateral vector from the center of the autocorrelation kernel
            x_vector = x_vector/max(real(x_vector)); % WE CHOOSE MAX AMZ
%             x_vector = x_vector/real(x_vector(window(2)));
            %center_x = x_vector(window(2));
            x_lag = x_vector(window(2) - d_n);
            %k_x(ii,jj) = sqrt( c_x*( 1 - real(x_lag)/real(center_x) ) ); % Lateral wavenumber matrix calculation with the approximation
            k_x(ii,jj) = sqrt( c_x*( 1 - real(x_lag) ) ); % Lateral wavenumber matrix calculation with the approximation

        end
    end

    k = (k_z + k_x)/2;
    sws_matrix = (2*pi*f_v)./k;
    sws_matrix = real(sws_matrix);
end

% function sws_matrix = sws_generator(u, window, f_v, d_n, dinf, og_size, est_z, est_x)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Function that yields the shear wave speed of a region with the 
% % approximation method. 
% % 
% % Inputs:  u           - 2D region of interest to evaluate
% %          window      - vector with the size of the window kernel
% %          f_v         - vibration freqeuncy %700 
% %          d_n         - Number of pixels to define the lag
% %          dinf        - structure that contains the spatial resolutions
% %          og_size     - vector containing the original size of the data
% %          est_z       - Axial estimator (usually 10)
% %          est_x       - Lateral estimator (usually 5)
% %
% % Outputs: sws_matrix  - Resulting shear wave speed matrix
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% 
%     res_z = dinf.dz; % Axial resolution
%     res_x = dinf.dx; % Lateral resolution
%     delta_z = d_n * res_z; % Axial lag
%     delta_x = d_n * res_x; % Lateral lag
%     c_z = est_z / (delta_z^2); % Axial constant 
%     c_x = est_x / (delta_x^2); % Lateral constant 
%     T_A = xcorr2(ones(window)); % Correction triangle
% 
%     % Initializations of the wavenumber and SWS matrices
%     k_z = zeros(length(1:2:og_size(1)), length(1:2:og_size(2)));
%     k_x = k_z;
%     sws_matrix = k_x;
%     step = 2;
% 
%     %% Begin
% 
%     search_z = 1:step:og_size(1) - window(1) + 1; % Ensure the window stays within bounds
%     for ii = 1:length(search_z)
% 
%         search_x = 1:step:og_size(2) - window(2) + 1; % Ensure the window stays within bounds
%         parfor jj = 1:length(search_x)
% 
%             % Window kernel 
%             area = u(search_z(ii):search_z(ii) + window(1) - 1, ...
%                      search_x(jj):search_x(jj) + window(2) - 1); 
% 
%             auto_matrix = xcorr2(area);
%             auto_matrix = auto_matrix ./ T_A; % Discrete function correction
% 
%             % Axial calculations
%             z_vector = auto_matrix(:, window(2))'; % Axial vector from the center of the autocorrelation kernel
%             z_vector = z_vector / max(real(z_vector)); % Normalize to max value
%             z_lag = z_vector(window(1) - d_n); % Calculate lag
%             k_z(ii, jj) = sqrt(c_z * (1 - real(z_lag))); % Axial wavenumber calculation
% 
%             % Lateral calculations
%             x_vector = auto_matrix(window(1), :); % Lateral vector from the center of the autocorrelation kernel
%             x_vector = x_vector / max(real(x_vector)); % Normalize to max value
%             x_lag = x_vector(window(2) - d_n); % Calculate lag
%             k_x(ii, jj) = sqrt(c_x * (1 - real(x_lag))); % Lateral wavenumber calculation
% 
%         end
%     end
% 
%     % Final shear wave speed matrix calculation
%     k = (k_z + k_x) / 2;
%     sws_matrix = (2 * pi * f_v) ./ k;
%     sws_matrix = real(sws_matrix);
% end
