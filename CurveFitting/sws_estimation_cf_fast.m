function [Kx, Kz, Rx, Rz, K1d, R1d] = sws_estimation_cf_fast(vz, win, dx, dz, correc, og_size)
% function [Kx, Kz, Rx, Rz] = sws_estimation_cf(vz, win, dx, dz, correc)
 % Combined function to estimate wavenumber and R-squared values
    % across both axial and lateral directions using parallel processing.
    
    % if nargin < 5
    %     correc = xcorr2(ones(win(1), win(2))); % Default correction if not provided
    % end

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
    

    step = 2;
    Rx = zeros(length(1:step:og_size(1)),length(1:step:og_size(2)));
    Rz = zeros(length(1:step:og_size(1)),length(1:step:og_size(2)));
    Kx = zeros(length(1:step:og_size(1)),length(1:step:og_size(2)));
    Kz = zeros(length(1:step:og_size(1)),length(1:step:og_size(2)));

    R1d = zeros(length(1:step:og_size(1)),length(1:step:og_size(2)));
    K1d = zeros(length(1:step:og_size(1)),length(1:step:og_size(2)));


    % extra save axis for faster computation 
    % EMZ style 

    % Axis of the autocorrelation kernels
    % z_axis = linspace(-(win(1)-1),win(1)-1,size(correc,1))*res_z;
    % x_axis = linspace(-(win(2)-1),win(2)-1,size(correc,2))*res_x;
    y_v = (-floor(size(correc,1)/2):1:floor(size(correc,1)/2))*dz;
    x_v = (-floor(size(correc,2)/2):1:floor(size(correc,2)/2))*dx;

    
    search_z = 1:step:og_size(1); % To save time, the windowed matrix is evaluated in points separated by a step of 2
    search_x = 1:step:og_size(2);


    
    % Double loop across both axial and lateral directions
    for i = 1:length(search_z)
        search_x = 1:step:og_size(2);
        % for j = 1:length(ev_al_pos_x)
        parfor j = 1:length(search_x)
            % Extract relevant submatrix from vz for current position
            % sub_vz = vz(ev_al_pos_z(i) + search_area_z, ev_al_pos_x(j) + search_area_x);
            
            sub_vz = vz(search_z(i):search_z(i)+win(1)-1,search_x(j):search_x(j)+win(2)-1); % Window kernel EMZ style
            
            % Call local loop function for wavenumber estimation
            [Raxial, Rlateral, k_axial, k_lateral, k_1d, R_1d] = test_reverb_cf_fast(sub_vz,dx,dz,mode,correc, y_v, x_v);

            % Store results
            Rx(i, j) = Rlateral;
            Rz(i, j) = Raxial;
            Kx(i, j) = k_lateral;
            Kz(i, j) = k_axial;
            R1d(i,j) = R_1d;
            K1d(i,j) = k_1d;
        end
    end
end

