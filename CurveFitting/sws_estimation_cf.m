function [Kx, Kz, Rx, Rz] = sws_estimation_cf(vz, win, dx, dz, correc)
% function [Kx, Kz, Rx, Rz] = sws_estimation_cf(vz, win, dx, dz, correc)
 % Combined function to estimate wavenumber and R-squared values
    % across both axial and lateral directions using parallel processing.
    
    if nargin < 5
        correc = xcorr(ones(win(1), win(2))); % Default correction if not provided
    end

    mode = 1; %1 for xcorr, 2 for w-k theorem
    
    [m, n] = size(vz);
    M = win(1);
    N = win(2);
    
    % Allow positions in both axial and lateral directions
    al_pos_z = round(M / 2):m - round(M / 2); % Axial positions
     
    ev_index_z = 1:floor(length(al_pos_z)); % even index vector 
    ev_al_pos_z = al_pos_z(ev_index_z);  % even positions from allow locations in axial direction
    search_area_z = -round(M / 2) + 1:round(M / 2) - 1; % Axial search area

    al_pos_x = round(N / 2):n - round(N / 2); % Lateral positions
    ev_index_x = 1:floor(length(al_pos_x)); % even index vector 
    ev_al_pos_x = al_pos_x(ev_index_x); % even positions from allow locations in lateral direction

    search_area_x = -round(N / 2) + 1:round(N / 2) - 1; % Lateral search area
    
    Rx = zeros(length(al_pos_z), length(al_pos_x)); % Preallocate results
    Rz = zeros(length(al_pos_z), length(al_pos_x));
    Kx = zeros(length(al_pos_z), length(al_pos_x));
    Kz = zeros(length(al_pos_z), length(al_pos_x));
    
    % Double loop across both axial and lateral directions
    for i = 1:length(ev_al_pos_z)
        parfor j = 1:length(ev_al_pos_x)
            % Extract relevant submatrix from vz for current position
            sub_vz = vz(ev_al_pos_z(i) + search_area_z, ev_al_pos_x(j) + search_area_x);
            
            % Call local loop function for wavenumber estimation
            [Raxial, Rlateral, k_axial, k_lateral] = test_reverb_cf(sub_vz,dx,dz,mode,correc);

            % Store results
            Rx(i, j) = Rlateral;
            Rz(i, j) = Raxial;
            Kx(i, j) = k_lateral;
            Kz(i, j) = k_axial;
        end
    end
end

