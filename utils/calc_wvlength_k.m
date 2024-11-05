function [wavelength, k] = calc_wvlength_k(sws, freq, dinf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [wavelength, k] = calc_wvlength_k(sws, freq, dinf)
% Function to calculate the wavelength in pixels and the wave number
%
% Inputs:
%   sws - Shear wave speed (in m/s)
%   freq - Frequency of the wave (in Hz)
%   dinf.dz - Axial resolution (in meters per pixel)
%   dinf.dx - Lateral resolution (in meters per pixel)
%
% Outputs:
%   wavelength - Wavelength in pixels along the axial direction
%       wavelength.meters  [m]
%       wavelength.pix_axi [pixels]
%       wavelength.pix_lat [pixels]
%   k - Wave number        [rad/m]
% AUTHORs: Edmundo Arom Miranda 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate wavelength in meters
    wavelength.meters = sws / freq;
    
    % Convert wavelength to pixels
    wavelength.pix_axi = round( wavelength.meters / dinf.dz );
    wavelength.pix_lat = round( wavelength.meters / dinf.dx );

    % Ensure wavelength_ax is an odd number
    if mod(wavelength.pix_axi, 2) == 0
        wavelength.pix_axi = wavelength.pix_axi + 1; % Make it odd by adding 1 if it's even
    end
    if mod(wavelength.pix_lat, 2) == 0
        wavelength.pix_lat = wavelength.pix_lat + 1; % Make it odd by adding 1 if it's even
    end
    
    % Calculate wave number (in radians per meter)
    k = 2 * pi * freq / sws;

    % Print results for visualization
    % fprintf(['SWS: %.2f [m/s] | Freq: %d [Hz] | K: %.2f [rad/m] | wvlength:  %.2f [mm] ' ...
    %     '| Window size Ax: %d | La: %d \n'], ...
    %  sws, freq, round(k), round(1e3*wavelength.meters), round(wavelength.pix_axi), round(wavelength.pix_lat)  );
    % 
    fprintf(['SWS: %.2f [m/s] | Freq: %d [Hz] | K: %.2f [rad/m] = %.4f [1/mm] | wvlength:  %.2f [mm] ' ...
        '| Window size Ax: %d | La: %d \n'], ...
     sws, freq, k, 1e-3*k/(2*pi),round(1e3*wavelength.meters), round(wavelength.pix_axi), round(wavelength.pix_lat)  );

end
