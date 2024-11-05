function u_Filtered = spatial_filtering(u_complex,dinf,f_v,cs_min,cs_max)
% function u_Filtered = spatial_filtering(u_complex,dinf,f_v,cs_min,cs_max)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %A 2-D FIR Bandpass Filter in the spatial domain with cutoff frequencies
% %set by cs_min and cs_max according to the relation k = 2*pi*f/c
% % 
% % Inputs:  u_complex   - Complex 2D matrix with temporal information
% %          f_v         - Vibration freqeuncy
% %          dinf        - Structure that contains the spatial resolutions
% %          cs_min      - Minimum wave velocity
% %          cs_max      - Maximim wave velocity
% % Outputs: u_Filtered  - Filtered data. Ready to apply the SWS estimator
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

    %%
    kl = (2*pi*f_v/cs_max); kh = (2*pi*f_v/cs_min);
    title_str = sprintf('2D Filt c_{H}=%.2f, k_{l}=%.2f rad/m | c_{l}=%.2f, k_{H}=%.2f rad/m', cs_max, kl, cs_min, kh);
    sigma = 300; % Filter order

    Fs_x = 1/dinf.dx;  % Sampling spatial frequency [px/m]
    Fs_z = 1/dinf.dz;  % Sampling spatial frequency [px/m]
    [KX, KZ] = freqspace(size(u_complex),'meshgrid');

    KX_m = KX*(Fs_x/2); % [1/m]
    KZ_m = KZ*(Fs_z/2); % [1/m]

    KX_rad = KX_m*2*pi; % [rad/m]
    KZ_rad = KZ_m*2*pi; % [rad/m]
    
    Hd = ones(size(u_complex)); 
    r = sqrt(KX_rad.^2 + KZ_rad.^2);
    Hd((r<kl*2)|(r>kh*2)) = 0;
    win = fspecial('gaussian',size(u_complex),sigma); 
    win = win ./ max(win(:));  % Make the maximum window value be 1.
    h = fwind2(Hd,win);        

    u_Filtered = filter2(h,u_complex); % Filtered signal

  %% PLOT FILTER EMZ VERSION
    kx_rad = KX_rad(1,:); %vector rad/m
    kz_rad = KZ_rad(:,1); %vector rad/m

    kx_m = KX_m(1,:); %vector 1/m
    kz_m = KZ_m(:,1); %vector 1/m
 
    H = freqz2(h, size(h));
    % H = fftshift( fft2(h) );

    figure,
    set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.55])
    tiledlayout(1,2)
    sgtitle(title_str);

    h1 = nexttile;
    imagesc(kx_rad, kz_rad, abs(H)), colorbar, colormap("parula"), grid;
    axis("image")
    xlabel('Kx (rad/m)'), ylabel('Ky (rad/m)')
    title('|H|')
    xlim([-5000 5000])
    ylim([-5000 5000])
    xticks([-5000 -2500 0 2500 5000])
    yticks([-5000 -2500 0 2500 5000])

    h1 = nexttile;
    imagesc(kx_m*1e-3, kz_m*1e-3, abs(H)), colorbar, colormap("parula"), grid;
    axis("image")
    xlabel('Kx (1/mm)'), ylabel('Ky (1/mm)')
    title('|H|')
    xlim([-2 2])
    ylim([-2 2])

% 
%     %% Plot section
%     % Plot filter
%     window = [100 100];
%     h_f = fftshift(fft2(h));
%     offset_z = floor(window(1)/2); offset_x = floor(window(2)/2);
%     center_z = floor(size(u_complex,1)/2);center_x = floor(size(u_complex,2)/2);
%     h_p = h_f(center_z-offset_z:center_z+offset_z-1,center_x-offset_x:center_x+offset_x-1);
%     z_axis = linspace(-offset_z,offset_z,window(1))*(Fsz/size_vector(1))*1e-3;
%     x_axis = linspace(-offset_x,offset_x,window(2))*(Fsx/size_vector(2))*1e-3;
% 
%     BPFilter = figure;
%     imagesc(x_axis,z_axis,abs(h_p)), colormap jet
%     xlabel('Spatial freq. [1/mm]'), ylabel('Spatial freq. [1/mm]')
%     title('2-D FIR Bandpass filter')
%     set(gca,'Fontname','Times')
%     % saveas(BPFilter,'FIR-2DBPF.png')
% 
%     % Plot input signal
%     plot_u = fftshift(fft2(u_complex));
%     plot_u = plot_u(center_z-offset_z:center_z+offset_z-1,center_x-offset_x:center_x+offset_x-1);
%     SFSIS = figure;
%     imagesc(x_axis,z_axis,abs(plot_u)), colormap jet
%     xlabel('Spatial freq. [1/mm]'), ylabel('Spatial freq. [1/mm]')
%     title('Spatial Frequency Spectrum of input signal')
%     set(gca,'Fontname','Times')
%     % saveas(SFSIS,'SpatFreqSpecInput.png')
% 
%     % Ploting Filter & Signal Profile Windows
%     kernel_size = 60;
%     v_ax = fftshift(fft(u_complex(:,center_z)'));
%     v_ax=v_ax./max(v_ax(:));
%     v_ax = v_ax(center_z-floor(kernel_size/2):center_z+floor(kernel_size/2));
%     v_lat = fftshift(fft(u_complex(:,center_x)));
%     v_lat=v_lat./max(v_lat(:));
%     v_lat = v_lat(center_x-floor(kernel_size/2):center_x+floor(kernel_size/2));
%     f_ax = fftshift(fft(h(:,center_z)'));
%     f_ax=f_ax./max(f_ax(:));
%     f_ax = f_ax(center_z-floor(kernel_size/2):center_z+floor(kernel_size/2));
%     f_lat = fftshift(fft(h(:,center_x)));
%     f_lat=f_lat./max(f_lat(:));
%     f_lat = f_lat(center_x-floor(kernel_size/2):center_x+floor(kernel_size/2));
% 
%     z_a = (-floor(kernel_size/2):floor(kernel_size/2))*(Fsz/size_vector(1))*1e-3;
%     x_a = (-floor(kernel_size/2):floor(kernel_size/2))*(Fsx/size_vector(1))*1e-3;
% 
%     filt_data_ax = figure;
%     plot(z_a,abs(v_ax),'k','Linewidth',1), hold on
%     plot(z_a,abs(f_ax),'k','Linewidth',2), hold off;
%     xlabel('Spatial Freq.[1/mm]'), ylabel('Normalized Magnitude')
%     title('Axial Profile'),legend('Data','Filter')
%     set(gca,'Fontname','Times')
%     % saveas(filt_data_ax,'AxProf.png')
% 
%     filt_data_lat = figure;
%     % plot(x_a,abs(v_lat),'k','Linewidth',1,z_a,abs(f_lat),'k','Linewidth',2)
%     plot(x_a,abs(v_lat),'k','Linewidth',1), hold on;
%     plot(z_a,abs(f_lat),'k','Linewidth',2), hold off;
%     xlabel('Spatial Freq.[1/mm]'), ylabel('Normalized Magnitude')
%     title('Lateral Profile'),legend('Data','Filter')
%     set(gca,'Fontname','Times')
%     % saveas(filt_data_lat,'LatProf.png')
% end