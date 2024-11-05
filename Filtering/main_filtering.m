%%
clc, clear all, close all
addpath(genpath(pwd))

sws2krad = @(sws, freq) 2*pi*freq / sws; 
krad2sws = @(krad, freq) 2*pi*freq / krad;

sws2kmm =  @(sws, freq) 1e-3*freq / sws; 
kmm2sws = @(kmm, freq) freq / (kmm*1e3);
%%
% load('D:\emirandaz\rswe\data_lim\12-09\300_360_400\12conc\RSWE_12conc_0500Hz.mat')
% load('D:\emirandaz\rswe\data_lim\16-09\300_360_400\RSWE_4conc_0500Hz.mat')

shaker_data = true;
if shaker_data
load('D:\emirandaz\rswe\data_lim\shakers\shakers-500.mat') % data shakers? 
end

freq = 500;
dinf.samplesRF = 2048; 
dinf.PRFe = 10*300; % dopPRF (LIM DATA) or dinf.PRFe

dinf.fc = Trans.frequency*1e6; % optional
dinf.c0 = 1540; % optional
dinf.wl = dinf.c0/dinf.fc; % optional
dinf.dz = PData(1).PDelta(3)*dinf.wl; % optional
dinf.dx = PData(1).PDelta(1)*dinf.wl; % optional

dinf.Tprf = 1/dinf.PRFe; 
dinf.Troundtrip = (dinf.samplesRF/2)*(1/dinf.fc);

dinf.maxDepth = dinf.Troundtrip*dinf.c0/2;
dinf.fs = dinf.PRFe; 
dinf.offset_y = 0; 
dinf.offset_z = 0;
dinf.num_angles = 1;


% shaker data 
if (shaker_data) 
    IQFrame = IQData(:,:,end); 
    IQData1 = IQData;
else 
    IQFrame = IQData1(:,:,end); % last frame
end
Bmode = 20*log10(abs(IQFrame(:,:,1)));
Bmode = Bmode - max(Bmode(:));

xdim = linspace(-dinf.dx*size(Bmode,2)/2,dinf.dx*size(Bmode,2)/2,size(Bmode,2)); 
% xdim = linspace(-dinf.Bdx*size(Bmode,2)/2,dinf.Bdx*size(Bmode,2)/2,size(Bmode,2)); % data2 LIM
ydim = linspace(0,dinf.dz*size(Bmode,1),size(Bmode,1));
x = xdim;
z = ydim;

[u,dinf] = pv_cal(IQData1,dinf, dinf.num_angles); 
u_new = signal_period(freq, dinf.PRFe, u);  

%% FGD EMZ ENHANCED
range_dB = [-120 0];
cs_min = 1; cs_max = 8; 
Fs = dinf.PRFe;
u_complex = temporal_extraction(u_new,freq,Fs);
u_Filtered = spatial_filtering(u_complex,dinf,freq,cs_min,cs_max);
%
% FFT2 plots

[KX, KZ] = freqspace(size(u_complex),'meshgrid');
Fs_x = 1/dinf.dx;
Fs_z = 1/dinf.dz;
KX_m = KX*(Fs_x/2); % [1/m]
KZ_m = KZ*(Fs_z/2); % [1/m]
KX_rad = KX_m*2*pi; % [rad/m]
KZ_rad = KZ_m*2*pi; % [rad/m]

kx_m = KX_m(1,:); 
kz_m = KZ_m(:,1);
kx_rad = KX_rad(1,:);
kz_rad = KZ_rad(:,1);

space_map = u_complex; 
kmap = fftshift( fft2(space_map) );
abs_kmap = abs(kmap);
abs_kmap_normalized = abs_kmap / max(abs_kmap(:));
abs_kmap_dB = 20 * log10(abs_kmap_normalized);

figure,
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.85 0.75])
tiledlayout(2,3)
sgtitle('')

nexttile
imagesc( x*1e3, z*1e3, angle(space_map)  );
axis('image')
title('u complex')
colorbar, 
% title(['k = ', num2str(wvnumber)])
xlabel('Lateral (mm)'), ylabel('Depth (mm)')

nexttile
imagesc(kx_rad, kz_rad, abs_kmap_dB, range_dB), colorbar, colormap("parula")
axis('image')
% imagesc(kx, ky, Hd), colorbar
xlabel('Kx (rad/m)'), ylabel('Ky (rad/m)')
% title('|H|')
a=colorbar;
a.Label.String = 'Norm. Power (dB)';
xlim([-5000 5000])
ylim([-5000 5000])

nexttile
imagesc(kx_m*1e-3, kz_m*1e-3,abs_kmap_dB, range_dB), colorbar, colormap("parula")
axis('image')
xlabel('Kx (1/mm)'), ylabel('Ky (1/mm)')
% title('|H|')
a=colorbar;
a.Label.String = 'Norm. Power (dB)';
xlim([-2 2])
ylim([-2 2])


space_map = u_Filtered; 
% Perform the 2D FFT and shift the zero-frequency component to the center
kmap = fftshift( fft2(space_map) );
abs_kmap = abs(kmap);
abs_kmap_normalized = abs_kmap / max(abs_kmap(:));
abs_kmap_dB = 20 * log10(abs_kmap_normalized);

nexttile
imagesc( x*1e3, z*1e3, angle(space_map)  );
title('u filtered')
axis('image')
colorbar, 
% title(['k = ', num2str(wvnumber)])
xlabel('Lateral (mm)'), ylabel('Depth (mm)')
% xlim([-5000 5000])
% ylim([-5000 5000])

nexttile
imagesc(kx_rad, kz_rad, abs_kmap_dB, range_dB), colorbar, colormap("parula")
axis('image')
% imagesc(kx, ky, Hd), colorbar
xlabel('Kx (rad/m)'), ylabel('Ky (rad/m)')
% title('|H|')
a=colorbar;
a.Label.String = 'Norm. Power (dB)';
xlim([-5000 5000])
ylim([-5000 5000])

nexttile
imagesc(kx_m*1e-3, kz_m*1e-3,abs_kmap_dB, range_dB), colormap("parula")
axis('image')
xlabel('Kx (1/mm)'), ylabel('Ky (1/mm)')
% title('|H|')
a=colorbar;
a.Label.String = 'Norm. Power (dB)';
xlim([-2 2])
ylim([-2 2])

%% JUVENAL EMZ ENHANCED
[M_filtered, M_complex] = fun_JO_v2(u_new, freq, dinf, cs_min, cs_max);

space_map = M_complex; 
kmap = fftshift( fft2(space_map) );
abs_kmap = abs(kmap);
abs_kmap_normalized = abs_kmap / max(abs_kmap(:));
abs_kmap_dB = 20 * log10(abs_kmap_normalized);

figure,
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.85 0.75])
tiledlayout(2,3)
sgtitle('')

nexttile
imagesc( x*1e3, z*1e3, angle(space_map)  );
axis('image')
title('u complex')
colorbar, 
% title(['k = ', num2str(wvnumber)])
xlabel('Lateral (mm)'), ylabel('Depth (mm)')

nexttile
imagesc(kx_rad, kz_rad, abs_kmap_dB, range_dB), colorbar, colormap("parula")
axis('image')
% imagesc(kx, ky, Hd), colorbar
xlabel('Kx (rad/m)'), ylabel('Ky (rad/m)')
% title('|H|')
a=colorbar;
a.Label.String = 'Norm. Power (dB)';
xlim([-5000 5000])
ylim([-5000 5000])

nexttile
imagesc(kx_m*1e-3, kz_m*1e-3,abs_kmap_dB, range_dB), colorbar, colormap("parula")
axis('image')
xlabel('Kx (1/mm)'), ylabel('Ky (1/mm)')
% title('|H|')
a=colorbar;
a.Label.String = 'Norm. Power (dB)';
xlim([-2 2])
ylim([-2 2])


space_map = M_filtered; 
% Perform the 2D FFT and shift the zero-frequency component to the center
kmap = fftshift( fft2(space_map) );
abs_kmap = abs(kmap);
abs_kmap_normalized = abs_kmap / max(abs_kmap(:));
abs_kmap_dB = 20 * log10(abs_kmap_normalized);

nexttile
imagesc( x*1e3, z*1e3, angle(space_map)  );
title('u filtered')
axis('image')
colorbar, 
% title(['k = ', num2str(wvnumber)])
xlabel('Lateral (mm)'), ylabel('Depth (mm)')
% xlim([-5000 5000])
% ylim([-5000 5000])

nexttile
imagesc(kx_rad, kz_rad, abs_kmap_dB, range_dB), colorbar, colormap("parula")
axis('image')
% imagesc(kx, ky, Hd), colorbar
xlabel('Kx (rad/m)'), ylabel('Ky (rad/m)')
% title('|H|')
a=colorbar;
a.Label.String = 'Norm. Power (dB)';
xlim([-5000 5000])
ylim([-5000 5000])

nexttile
imagesc(kx_m*1e-3, kz_m*1e-3,abs_kmap_dB, range_dB), colormap("parula")
axis('image')
xlabel('Kx (1/mm)'), ylabel('Ky (1/mm)')
% title('|H|')
a=colorbar;
a.Label.String = 'Norm. Power (dB)';
xlim([-2 2])
ylim([-2 2])

%
%% LIM ENHANCED
f_tol = 50;
[u_new, Frames0, Frames1] = u_filt(u, freq, f_tol, dinf, cs_min, cs_max);  

%% 

%
% 
% name_path = 'C:\Users\u_imagenes\Desktop\Fernando';
% 
% for aa = 1:18 % limit term must be varied to the number of files in the folder
%     path1=[name_path '\File_' num2str(mm)];
%     cd(path1)
%     load('pv_cal.mat')
%     for f = freq
%         mkdir([path1,'\data_freq_',num2str(f)])                                  
%         cd([path1,'\data_freq_',num2str(f)])
%         cs_min = 0.7; cs_max = 5; 
%         Fs = dinf.PRFe;
%         u_complex = temporal_extraction(u,f,F_s);
%         u_Filtered = spatial_filtering(u_complex,dinf,f_v,cs_min,cs_max);
%         save('data_filt.mat',u_filtered);
%     end
% end


%%
