% Read data k-wave generated by cluster
% Verify reverberation by kmap analys
% Updated November

clear all, clc;

%% Some utils   

sws2krad = @(sws, freq) 2*pi*freq / sws; 
krad2sws = @(krad, freq) 2*pi*freq / krad;

sws2kmm =  @(sws, freq) 1e-3*freq / sws; 
kmm2sws = @(kmm, freq) freq / (kmm*1e3);

%% (i) 3 attempt
%% LOAD DATA
pathdata = 'P:\emiranda\proj\rswe\kwave_v1\simRSWF\';
% namedata = 'RSWF_source3_homo'; %grid size 50x50x50

% good simulations 60x60x60
% namedata = 'RSWF_homo_s3_f1000';
namedata = 'RSWF_homo_s12_f1000';
% namedata = 'RSWF_homo_s50_f1000';
% namedata = 'RSWF_homo_s200_f1000';

mykwave = load(fullfile(pathdata, namedata));

%% EXTRACT ARRAY PART VELOCITY AND AXIS
array2D = squeeze(mykwave.uz_s_array(10,:,:,:));

xx = mykwave.kgrid.x_vec;
zz = mykwave.kgrid.z_vec;

%% VIDEO CHECK
% 
% max_uz = max(mykwave.uz_s_array(:));
% min_uz = min(mykwave.uz_s_array(:));
% 
% figure,
% for ii = 2000:5:3000
%     %im=squeeze((uz_s(10,:,:,ii)-min_uz)/(max_uz-min_uz))';
%     im = squeeze((mykwave.uz_s_array(10,:,:,ii)))';
%     imagesc(xx*1E3, zz*1E3, im);
%     xlabel('x (mm)'); ylabel('z (mm)');
%     axis image;
%     % clim([-0.01 0.01]); %modify limits to see waves
%     clim([-10E-5 10E-5]); %modify limits to see waves
%     title(['SWS - XZ plane - Frame ' num2str(ii)]);
%     grid on;
%     colormap('jet')
%     colorbar;
%     drawnow
%     pause(0.001)
% 
% end

%% extract pixel
[m, n, p, r] = size(mykwave.uz_s_array);
pixel = mykwave.uz_s_array(ceil(m/2), ceil(n/2), ceil(p/2), :);
pixel = squeeze(pixel);

t_vec = (0 : mykwave.kgrid.Nt - 1)*mykwave.kgrid.dt;
% figure, 
% plot(t_vec*1E3, pixel), grid on;
% xlabel('time (ms)')


%% Extract steadystate

uz_steadyState = array2D(:,:,2000:end); % last 1000 frames

[m,n,p] = size(uz_steadyState);


%% FFT TEMPORAL EXTRACTION 

sws_gt = 1;  % [m/s]
freq = 1000; % [Hz]
Fs = 1/mykwave.kgrid.dt;

u_complex = temporal_extraction(uz_steadyState, freq, Fs);

% 
% figure, 
% imagesc(xx, zz, angle(u_complex))

%% K-SPACE PLOT NORMALIZED
fontSize = 20;

[KX, KZ] = freqspace(size(u_complex),'meshgrid');

Fs_x = 1/mykwave.kgrid.dx;
Fs_z = 1/mykwave.kgrid.dy;

KX_m = KX*(Fs_x/2); % [1/m]
KZ_m = KZ*(Fs_z/2); % [1/m]
KX_rad = KX_m*2*pi; % [rad/m]
KZ_rad = KZ_m*2*pi; % [rad/m]
    
kx_m = KX_m(1,:); 
kz_m = KZ_m(:,1);
kx_rad = KX_rad(1,:);
kz_rad = KZ_rad(:,1);

space_map = angle( u_complex ); 
kmap = fftshift( fft2(space_map) );

abs_kmap = abs(kmap);
abs_kmap_normalized = abs_kmap / max(abs_kmap(:));
abs_kmap_dB = 20 * log10(abs_kmap_normalized);
range_dB = [-60 0];

%%%%%%%%%%%%%%% IDEAL K VALUES %%%%%%%%%%%%%%%
k_gt_rad = sws2krad(sws_gt, freq);
k_gt_mm  = sws2kmm(sws_gt, freq);

% To create a circle plot
theta = linspace(0, 2*pi, 100);
x_circle_rad = k_gt_rad * cos(theta);
y_circle_rad = k_gt_rad * sin(theta);

x_circle_mm = k_gt_mm * cos(theta);
y_circle_mm = k_gt_mm * sin(theta);
%%%%%%%%%%%%%%% IDEAL K VALUES %%%%%%%%%%%%%%%

figure,
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.85 0.75])
tiledlayout(1,3)
sgtitle('')

nexttile
imagesc( xx*1e3, zz*1e3, space_map  );
axis('image')
title('u complex')
colorbar, 
% title(['k = ', num2str(wvnumber)])
xlabel('Lateral (mm)'), ylabel('Depth (mm)')

nexttile
imagesc(kx_rad, kz_rad, abs_kmap_dB, range_dB), colorbar, colormap("parula")

hold on;
plot(x_circle_rad, y_circle_rad, 'r--', 'LineWidth', 2)  % Plot in white for visibility
hold off;

axis('image')
title('Kspace (rad/m)')
% imagesc(kx, ky, Hd), colorbar
xlabel('Kx (rad/m)'), ylabel('Ky (rad/m)')
% title('|H|')
a=colorbar;
a.Label.String = 'Norm. Power (dB)';
% xlim([-8000 8000])
% ylim([-8000 8000])

nexttile
imagesc(kx_m*1e-3, kz_m*1e-3,abs_kmap_dB, range_dB), colorbar, colormap("parula")

hold on;
plot(x_circle_mm, y_circle_mm, 'r--', 'LineWidth', 2)  % Plot in white for visibility
hold off;

axis('image')
title('Kspace (1/mm)')
xlabel('Kx (1/mm)'), ylabel('Ky (1/mm)')
% title('|H|')
a=colorbar;
a.Label.String = 'Norm. Power (dB)';
% xlim([-3 3])
% ylim([-3 3])

% Apply font size to all axes labels, titles, and colorbars in one line
set(findall(gcf,'Type','text'),'FontSize',fontSize); % For all text elements
set(findall(gcf,'Type','axes'),'FontSize',fontSize); % For axis labels and titles


%%
% % Initialize the result array for storing the FFT
% uz_steadyState_fft = zeros(m, n, p);
% 
% % Calculate the FFT along the third dimension (frame-by-frame) for each pixel
% for i = 1:m
%     for j = 1:n
%         % Perform FFT and FFT shift on the third dimension for each pixel
%         uz_steadyState_fft(i, j, :) = fftshift(fft(uz_steadyState(i, j, :)));
%     end
% end
% 
% %
% frame = uz_steadyState_fft(:,:,100);
% 
% figure, 
% imagesc(xx, zz, angle(frame))
% colormap("parula")

%%
figure,
histogram(abs_kmap_dB), 
title('Histogram Norm Kmap dB');

%%
%% K-SPACE PLOT
fontSize = 20;

[KX, KZ] = freqspace(size(u_complex),'meshgrid');

Fs_x = 1/mykwave.kgrid.dx;
Fs_z = 1/mykwave.kgrid.dy;

KX_m = KX*(Fs_x/2); % [1/m]
KZ_m = KZ*(Fs_z/2); % [1/m]
KX_rad = KX_m*2*pi; % [rad/m]
KZ_rad = KZ_m*2*pi; % [rad/m]
    
kx_m = KX_m(1,:); 
kz_m = KZ_m(:,1);
kx_rad = KX_rad(1,:);
kz_rad = KZ_rad(:,1);

space_map = angle( u_complex ); 
kmap = fftshift( fft2(space_map) );

abs_kmap = abs(kmap);
% abs_kmap_normalized = abs_kmap / max(abs_kmap(:));
abs_kmap_dB = 20 * log10(abs_kmap);
range_dB = [-0.8 65.8];

%%%%%%%%%%%%%%% IDEAL K VALUES %%%%%%%%%%%%%%%
k_gt_rad = sws2krad(sws_gt, freq);
k_gt_mm  = sws2kmm(sws_gt, freq);

% To create a circle plot
theta = linspace(0, 2*pi, 100);
x_circle_rad = k_gt_rad * cos(theta);
y_circle_rad = k_gt_rad * sin(theta);

x_circle_mm = k_gt_mm * cos(theta);
y_circle_mm = k_gt_mm * sin(theta);
%%%%%%%%%%%%%%% IDEAL K VALUES %%%%%%%%%%%%%%%

figure,
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.85 0.75])
tiledlayout(1,3)
sgtitle('')

nexttile
imagesc( xx*1e3, zz*1e3, space_map  );
axis('image')
title('u complex')
colorbar, 
% title(['k = ', num2str(wvnumber)])
xlabel('Lateral (mm)'), ylabel('Depth (mm)')

nexttile
imagesc(kx_rad, kz_rad, abs_kmap_dB, range_dB), colorbar, colormap("parula")
% imagesc(kx_rad, kz_rad, abs_kmap_dB), colorbar, colormap("parula")

hold on;
plot(x_circle_rad, y_circle_rad, 'r--', 'LineWidth', 2)  % Plot in white for visibility
hold off;

axis('image')
title('Kspace (rad/m)')
% imagesc(kx, ky, Hd), colorbar
xlabel('Kx (rad/m)'), ylabel('Ky (rad/m)')
% title('|H|')
a=colorbar;
a.Label.String = 'Power (dB)';
% xlim([-8000 8000])
% ylim([-8000 8000])

nexttile
imagesc(kx_m*1e-3, kz_m*1e-3,abs_kmap_dB, range_dB), colorbar, colormap("parula")
% imagesc(kx_m*1e-3, kz_m*1e-3,abs_kmap_dB), colorbar, colormap("parula")

hold on;
plot(x_circle_mm, y_circle_mm, 'r--', 'LineWidth', 2)  % Plot in white for visibility
hold off;

axis('image')
title('Kspace (1/mm)')
xlabel('Kx (1/mm)'), ylabel('Ky (1/mm)')
% title('|H|')
a=colorbar;
a.Label.String = 'Power (dB)';
% xlim([-3 3])
% ylim([-3 3])

% Apply font size to all axes labels, titles, and colorbars in one line
set(findall(gcf,'Type','text'),'FontSize',fontSize); % For all text elements
set(findall(gcf,'Type','axes'),'FontSize',fontSize); % For axis labels and titles
%%
figure,
histogram(abs_kmap_dB), 
title('Histogram Kmap dB');