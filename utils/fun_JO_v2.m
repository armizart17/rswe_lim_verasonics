function [M_filtered] = fun_JO_v2(u, freq, dinf)
% function [M_filtered] = fun_JO_v1(u, freq, dinf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple function for phase extraction based on Juvenal Ormachea code
% Inputs:  
%           u           : input particle velocity
%           freq        : Vibration freqeuncy
%           dinf        : Frequency range for the bandpass cuttoffs are +-
%                         2*f_band
% Outputs: 
%           M_filtered  : filtered particle velocity signal
% Author: Edited by EMZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % u = u(:,:,50:end);  % usually 50 first frames are corrupted
    freq = freq(1);
    Fs = dinf.PRFe;      % Sampling frequency

    [a, b, c] = size(u);
%     Ts = 1/Fs;                            % Sample time
    L = 2^nextpow2(c*3);                    % Length of signal
    df = Fs/L;                              % frequency resolution
    f1 = (-round(L/2):1:round(L/2)-1)*df;   % frequency axis
    %% Find location of frequencies
    cont = 1;
    ix = zeros(1,size(freq,2));
    for ii = freq
        [ ~, ix(cont) ] = min( abs( f1-ii ));
        cont = cont+1;
    end
    %% Magnitude and phase extraction  
    u2 = reshape(u,[a*b c]);
    M_complex = zeros(a*b,size(freq,2));
    
    for ii = 1: a*b
        F1 = fftshift(fft(u2(ii,:),L));
        M_complex(ii,:) = F1(ix);
    end
    M_complex = reshape(M_complex,[a b size(freq,2)]);
    
    sigma = 50; % CHANGE 300;
    %  Breast data
    cs_min = 0.045*(2*pi*freq).^0.25;
    cs_max = 0.33*(2*pi*freq).^0.35 /2; % I CHANGED THIS TO MAKE IT WORK
    % cs_max = 4.5;
    Fs1 = 1/dinf.dz;  % Sampling spatial frequency
    Fs2 = 1/dinf.dx;


    %%
[Ny, Nx, ~] = size(u);

% Define the spatial sampling intervals (assuming equal spacing)
dx = dinf.dx;  % Example spatial resolution in meters
dy = dinf.dz;       % Assume square grid for simplicity

% Compute the wavenumber ranges
kx = (-Nx/2:Nx/2-1) * (2*pi / (Nx * dx));
ky = (-Ny/2:Ny/2-1) * (2*pi / (Ny * dy));

% Create the wavenumber grid
[KX, KY] = meshgrid(kx, ky);




%%
    
    r = sqrt(KX.^2 + KY.^2);
    win = fspecial('gaussian',[a b],sigma);
    win = win ./ max(win(:));  % Make the maximum window value be 1.
    
    M_filtered = zeros(a,b,size(freq,2));
    cont = 1;
    for ii = freq
        Hd = ones([a b]);
        Hd1 = ones([a b]);
        kl = (2*pi*ii/cs_max(cont)); kh = (2*pi*ii/cs_min(cont));
        Hd((r<kl)|(r>kh)) = 0;
        figure, 
        imagesc( kx, ky, Hd  )
        xlabel('Kx (rad/m)'), ylabel('Ky (rad/m)')
        title('Hd'), colorbar;
        xlim([-20e3, 20e3])
        ylim([-20e3, 20e3])
        h = fwind2(Hd,win);  % Using the 2-D window, design the filter that best produces the desired frequency response
        M_filtered(:,:,cont) = filter2(h,M_complex(:,:,cont));
        cont = cont+1;
    end
    
    mask_ph = fspecial('average',[6 6]);
    M_filtered = imfilter(M_filtered, mask_ph,'replicate');
end

