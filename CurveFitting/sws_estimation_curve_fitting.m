function [Kx,Kz,Rx,Rz] = sws_estimation_curve_fitting(vz, win, dx , dz, correc)
% function [Kx,Kz,Rx,Rz] = sws_estimation_curve_fitting(vz, win, dx , dz, correc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local wavenumber estimator by Curve Fitting LIM version
% Inputs:  vz       - Complex matrix with magnitude and phase info from RSWField
%          window   - Size of the windows kernel.
%          dx       - Resolution of the x-axis (lateral) of the signal
%          dz       - Resolution of the z-axis (axial) of the signal
%          correct  - should be xcorr(ones(win(1), win(2))
% Outputs: Kx,Kz    - Spatial frequency or wavenumber (k)
%          Rx,Rz    - R - square (k)
% Author: E.A.M.Z. adapted from LIM code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If correc is not provided, calculate it using the default value
    if nargin < 5
        correc = xcorr(ones(win(1), win(2)));
    end
    
    [m,n] = size(vz);
    M = win(1); 
    N = win(2); 
        
    al_pos_z = round(M/2):m-round(M/2); % allow locations to use in axial direction
    ev_index_z = 1:floor(length(al_pos_z)); % even index vector 
    ev_al_pos_z = al_pos_z(ev_index_z);  % even positions from allow locations in axial direction
    search_area_z = -round(M/2)+1:round(M/2)-1; % search index from kernel:lateral
    
    % Preallocate matrices for Rx, Rz, Kx, Kz
    num_pos_z = length(ev_al_pos_z);
    Rx = zeros(num_pos_z, n);  % Lateral R-squared values
    Rz = zeros(num_pos_z, n);  % Axial R-squared values
    Kx = zeros(num_pos_z, n);  % Lateral wavenumbers
    Kz = zeros(num_pos_z, n);  % Axial wavenumbers

    for k = 1:length(ev_al_pos_z)
    %  2D cross-corralation of the reverberant particle velocity. 
    %  Bvv (from papers) are extracted at different directions (angles)
        [Raxial,Rlateral,K_axial,K_lateral] = sws_estimation_localloop_v1(vz(ev_al_pos_z(k)+search_area_z,:),dx,dz,N,n,correc);
        Rx(k,:) = Rlateral;
        Rz(k,:) = Raxial;
        Kx(k,:) = K_lateral;
        Kz(k,:) = K_axial;
    end 
end