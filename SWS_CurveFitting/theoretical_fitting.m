function [k_z,R_ax,k_x,R_lat,k,sws_matrix] = theoretical_fitting(u,window,f_v,dinf,og_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the shear wave speed of a region with the 
% curve fitting method. Its corresponding parameters of goodness of Fit
% (R^2) are also stored.
% 
% 
% Inputs:  u           - 2D region of interest to evaluate (windowed data)
%          window      - vector with the size of the window kernel
%          f_v         - vibration freqeuncy
%          dinf        - structure that contains the spatial resolutions
%          og_size     - vector containing the original size of the data
%
% Outputs: k_z         - Wavenumber matrix for the axial direction
%          R_ax        - Goodness of Fit matrix for axial direction
%          k_x         - Wavenumber matrix for the lateral direction
%          R_lat       - Goodness of Fit matrix for lateral direction
%          k           - Average wavenumber matrix 
%          sws_matrix  - Resulting shear wave speed matrix 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resolution
    T_A = xcorr2(ones(window)); % 2D Correction Triangle
    step = 2;
    %% Initializing vectors and matrixes
    k_z = zeros(length(1:step:og_size(1)),length(1:step:og_size(2))); 
    k_x = k_z; % matrixes containing the estimated wavenumber along each direction
    sws_matrix = k_x; % Matrix containing the calculated Shear Wave Speed
    R_ax = sws_matrix; 
    R_lat = R_ax;
        % Axis of the autocorrelation kernels
    z_axis = linspace(-(window(1)-1),window(1)-1,size(T_A,1))*res_z;
    x_axis = linspace(-(window(2)-1),window(2)-1,size(T_A,2))*res_x;
    
    %% Begin
    
    search_z = 1:step:og_size(1); % To save time, the windowed matrix is evaluated in points separated by a step of 2
    %
    search_x = 1:step:og_size(2);
    %
    for ii = 1:length(search_z)
        %%%
        search_x = 1:step:og_size(2);
        %%%
        parfor jj=1:length(search_x)
        % for jj=1:length(search_x)
        %for jj = 1:2:og_size(2)
            area = u(search_z(ii):search_z(ii)+window(1)-1,search_x(jj):search_x(jj)+window(2)-1); % Window kernel 
            auto_matrix = xcorr2(area); % Autocorrelation of the kernel
            auto_matrix = auto_matrix./T_A; % correction for discrete functions
            auto_matrix = auto_matrix./( auto_matrix( window(1),window(2) ) ); % Normalizing
            
            z_vector = real(auto_matrix(:,window(2))'); % Axial vector from the center of the autocorrelation kernel
            x_vector = real(auto_matrix(window(1),:)); % Lateral vector from the center of the autocorrelation kernel
            
            %% To exclude oscillations located far from the center (From Ormachea code - optional)
            track = 0.005;
%             [z_axis2,z_vector2,x_axis2,x_vector2] = reduce_axis(track,z_axis,z_vector,x_axis,x_vector);
            %%
            z_axis2 = z_axis;
            z_vector2 = z_vector;
            x_axis2 = x_axis;
            x_vector2 = x_vector;
            %% Curve fitting to the theorical Bessel Functions for both axial and lateral directions
            [result1, gof1] = rswe_besselfit_axial(z_axis2, z_vector2);
            [result2, gof2] = rswe_besselfit_lateral(x_axis2, x_vector2);
            
            % Storing the values of the wavenumbers and the GOF
            k_z(ii,jj) = result1.k;
            k_x(ii,jj) = result2.k;
            R_ax(ii,jj) = gof1(1).rsquare;
            R_lat(ii,jj) = gof2(1).rsquare;
            
        end
    end
    
    k = (k_z + k_x)/2;
    sws_matrix = (2*pi*f_v)./k;
    sws_matrix = real(sws_matrix);
end