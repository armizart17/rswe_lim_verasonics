% R-SWE Imaging: FAST TEST CF, MAOW, PG-Least Squares and PG-TV Estimators
% Description: 
% This code generates part veloc, and SWS by different estimators for R-SWE
% Saved in: Customizable
% Author: EMZ (based on LIM codes)
% myBmode = @(RF) 20*log10(abs(hilbert(RF))) - max(20*log10(abs(hilbert(RF(:)))));

% ASSUME TYPICAL DATA WILL BE SAVE LIKE THIS
% FolderName: 'dateofAcq' ------> i.e. '17-11'
% TypicalName:  'RSWE_idName.mat' ----->  i.e. 'RSWE_2lay_0500Hz.mat'
% idName = 2lay_0500Hz  General Recommended Structure: "typePhatom_freq(4dig)Hz"

% Description:
% First test for G. Flores and S. Romero data from UofR Feb 24th, 2025
% Data saved in "D:\emirandaz\rswe\data_lim\data_freq_300_Arom"

% SPECS
% SWS: 2.04 [m/s] | Freq: 300 [Hz] | K: 924.00 [rad/m] = 0.1471 [1/mm] | wvlength:  7.00 [mm] | Window size Ax: 41 | La: 21 
% SWS: 2.88 [m/s] | Freq: 300 [Hz] | K: 654.50 [rad/m] = 0.1042 [1/mm] | wvlength:  10.00 [mm] | Window size Ax: 57 | La: 29 

% Gilmer data

addpath(genpath(pwd));

%% Setup Directories and FolderName and Methods
clear all, clc
close all;

baseDir = 'D:\emirandaz\rswe\data_lim';     %%%% *CHANGE IT IF NEED IT %%%%
folderAcqName = 'data_freq_300_Arom' ;      %%%% *CHANGE IT IF NEED IT %%%%

dataDir = fullfile(baseDir, folderAcqName);

fprintf('==============================================================\n');
fprintf('Data Directory Path: %s\n', dataDir);
fprintf('==============================================================\n');

% outDir = ; % for Results 

methodCF    = false;
methodMAOW  = false;
methodPG_LS = true;
methodPG_TV = true;
%% Setup TypicalName of Data

% typName = input('Enter full name of data (i.e. RSWE_2lay_500Hz): ', 's');

% pos_ = strfind(typName, '_'); **
% posHz = strfind(typName, 'Hz'); **

% idName = typName(pos_(1)+1 : posHz+1); **
idName = 'v1';

% fprintf('==============================================================\n');
% fprintf('Data Full Name: %s\n', typName)
% fprintf('Data ID Name: %s\n', idName)
% fprintf('==============================================================\n');

% Extract the four digits before "Hz"
% freq = str2double( typName(posHz-4:posHz-1) ); % frequency usually represented by 4digits **
freq = 300;
%% LOAD DATA AND BMODE VALIDATION
% load(fullfile(dataDir, typName +".mat")); **

load(fullfile(dataDir,'pv_cal.mat' ))
load(fullfile(dataDir,'data_filt.mat' ))
load(fullfile(dataDir,'binaryMask_sws.mat'))

calc_wvlength_k(2.04, freq, dinf); %cs_min = 2.04m/s, cs_max = 2.88
calc_wvlength_k(2.88, freq, dinf); %cs_min = 2.04m/s, cs_max = 2.88

% SWS: 2.04 [m/s] | Freq: 300 [Hz] | K: 924.00 [rad/m] = 0.1471 [1/mm] | wvlength:  7.00 [mm] | Window size Ax: 41 | La: 21 
% SWS: 2.88 [m/s] | Freq: 300 [Hz] | K: 654.50 [rad/m] = 0.1042 [1/mm] | wvlength:  10.00 [mm] | Window size Ax: 57 | La: 29 

% dinf.fc = Trans.frequency*1e6;
% dinf.c0 = 1540; dinf.wl = dinf.c0/dinf.fc;
% dinf.dz = PData(1).PDelta(3)*dinf.wl;
% dinf.dx = PData(1).PDelta(1)*dinf.wl;
% dinf.samplesRF = 2048; dinf.PRFe = 10*500;
% dinf.Tprf = 1/dinf.PRFe; dinf.Troundtrip = (dinf.samplesRF/2)*(1/dinf.fc);
% dinf.maxDepth = dinf.Troundtrip*dinf.c0/2;
% dinf.fs = 5000; fs = 5000;dinf.offset_y = 0; 

% IQFrame = IQData1(:,:,end); % last frame
% Bmode = 20*log10(abs(IQFrame(:,:,1)));
% Bmode = Bmode - max(Bmode(:));

% dinf.num_angles = 1;
% dinf.offset_z = dinf.dz*size(IQFrame,1)/2;

xdim = linspace(-dinf.dx*size(Bmode,2)/2,dinf.dx*size(Bmode,2)/2,size(Bmode,2)); 
ydim = linspace(0,dinf.dz*size(Bmode,1),size(Bmode,1));

figure, 
imagesc(xdim*1e2, ydim*1e2, Bmode), 
colormap('gray'), 
xlabel('\bfLateral [cm]'),
ylabel('\bfAxial [cm]');
title(['B-mode ID: ', idName], 'Interpreter', 'none', 'FontWeight', 'bold')


%% FAST TEST KMAP
f_tol = 50;
cs_min = 1; cs_max = 5;
range_dB = [-120 0];
[u_new, Wave_z, Frames1, Frames0] = u_filt(u, freq, f_tol, dinf, cs_min, cs_max);  

figure, 
subplot(231), imagesc(abs(Wave_z)), title('Abs Wave_z')
subplot(234), imagesc(angle(Wave_z)), title('Angle Wave_z')

subplot(232), imagesc(abs(Frames1)), title('Abs Frames1')
subplot(235), imagesc(angle(Frames1)), title('Angle Frames1')

subplot(233), imagesc(abs(Frames0)), title('Abs Frames0')
subplot(236), imagesc(angle(Frames0)), title('Angle Frames0')

[KX, KZ] = freqspace(size(Wave_z),'meshgrid');
Fs_x = 1/dinf.dx; Fs_z = 1/dinf.dz;

KX_m = KX*(Fs_x/2); KZ_m = KZ*(Fs_z/2); % [1/m]
KX_rad = KX_m*2*pi; KZ_rad = KZ_m*2*pi; % [rad/m]

kx_m = KX_m(1,:); kz_m = KZ_m(:,1);
kx_rad = KX_rad(1,:); kz_rad = KZ_rad(:,1);

clear KX_m KZ_m KX_rad KZ_rad

space_map = Wave_z; 

kmap = fftshift( fft2(space_map) );
abs_kmap = abs(kmap);
abs_kmap_normalized = abs_kmap / max(abs_kmap(:));
abs_kmap_dB = 20 * log10(abs_kmap_normalized);

figure,
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.85 0.75])
tiledlayout(2,3)
sgtitle('')

nexttile
imagesc( xdim*1e3, ydim*1e3, angle(space_map)  );
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

%% PARTICLE VELOCITY ESTIMATION 

% [u,dinf] = pv_cal(IQData1,dinf, dinf.num_angles); 
% u = signal_period(freq, dinf.PRFe, u);   

% Temperal filtering process, a bandpass FIR filter is used around 
% +- 20 Hz the vibration frequency
cs_min = 0.7; % [m/s]
cs_max = 5;   % [m/s]
f_tol = 50;   % [Hz] tolerance for +/-2*f_tol
% This function uses spatial_fil_phase_extrac inside


% [u_new, Wave_z, Frames1] = u_filt(u, freq, f_tol, dinf, cs_min, cs_max);  
% u_new: video
% Wave_z: complex PV
% Frames1: onlyPhase PV



% CHECK PLOT Wave_z and FRAMES 1
% figure(), 
% sgtitle(['Wave_z ID: ', idName], 'Interpreter', 'none')
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
% subplot(121), 
% imagesc(xdim*1e2, ydim*1e2,abs(Wave_z)), title('Magn'), axis("tight");
% xlabel('Lateral [cm]'), ylabel('Axial [cm]')
% % set(gca, "clim", [1, 5])
% subplot(122), 
% imagesc(xdim*1e2, ydim*1e2, angle(Wave_z)), title('Phase'), axis("tight"), colormap('default')
% 
% figure(), 
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
% sgtitle(['Frames1 ID: ', idName], 'Interpreter', 'none')
% subplot(121), 
% imagesc(abs(Frames1)), title('Mang'), axis("tight");
% % set(gca, "clim", [1, 5])
% subplot(122), 
% imagesc(angle(Frames1)), title('Phase'), axis("tight"), colormap('default')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SWS ESTIMATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% VISUALIZATION SWS IMAGES DEFAULT

% Define default Visualization specs
visDefault.x = xdim;                % Lateral coordinates
visDefault.z = ydim;                % Axial coordinates
visDefault.BmodeFull = Bmode;       % B-mode image data
visDefault.caxis_bmode = [-60 0];   % Color axis limits for B-mode
visDefault.caxis_img = [0 4.5];     % Color axis limits for color (SWS)
visDefault.fact_transparency = 0.6; % Example transparency factor

%% CURVE FITTING (CF) "VERSION EMZ"
if methodCF

methodName = 'CF';
window = 101; 
w_kernel = [window, window];

% pv_field = Wave_z;  
pv_field = angle(exp_phase);
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

[k_z,R_ax,k_x,R_lat,k,sws_cf] = theoretical_fitting(mirror_frame,w_kernel,freq,dinf,og_size);

sws_cf_big = bigImg(sws_cf, Bmode); % resize to Bmode size

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizCF = Visualizer_v2(visDefault);

vizCF = vizCF.setROI(xdim, ydim, sws_cf_big);

titleName = strcat(methodName, ' ID-', idName);
vizCF = vizCF.setTitle(titleName);
vizCF = vizCF.setUnits('mm');

vizCF.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
end
%% CURVE FITTING sws_estimation_cf_fast
win = [71 71];
correc = xcorr2(ones(win(1), win(2))); % it is default 

pv_field = angle(exp_phase);
og_size = size(pv_field);
mirror_pv = padarray(pv_field, (win-1)/2, 'symmetric','both');

tic 
% [Kx,Kz,Rx,Rz] = sws_estimation_curve_fitting(mirror_pv, win, dx , dz, correc);
% [Kx,Kz,Rx,Rz] = sws_estimation_cf(mirror_pv, win, dx, dz, correc);

% [Kx,Kz,Rx,Rz, K1d, R1d] = sws_estimation_cf(mirror_pv, win, dx, dz, correc); % GILMER
[Kx,Kz,Rx,Rz,K1d,R1d] = sws_estimation_cf_fast(mirror_pv, win, dinf.dx, dinf.dz, correc, og_size); % EMZ
tt = toc;
fprintf('Time passed CF %.4f\n', tt)

%%
figure, 
sgtitle(['Win = ', num2str(win)])
wvnumber = 900;
axis_k = [0 10*wvnumber];
subplot(2,3,1), imagesc(real(Kx), axis_k), colorbar, title('Kx')
subplot(2,3,2), imagesc(real(Kz), axis_k), colorbar, title('Kz')
subplot(2,3,3), imagesc(real(K1d), axis_k), colorbar, title('K1d')

subplot(2,3,4), imagesc(real(Rx)), colorbar, title('Rx')
subplot(2,3,5), imagesc(real(Rz)), colorbar, title('Rz')
subplot(2,3,6), imagesc(real(R1d)), colorbar, title('R1d')

%
sws_cf_lat = 2*pi*freq./Kx;
sws_cf_axi = 2*pi*freq./Kz;
sws_cf_both = 2*pi*freq./( 0.5*(Kz+Kx)  );
sws_cf_1d = 2*pi*freq./K1d;

figure, 
sgtitle(['Win = ', num2str(win)])
subplot(2,2,1), imagesc(sws_cf_lat, visDefault.caxis_img), title('SWS lat'), colormap("turbo"), colorbar
subplot(2,2,2), imagesc(sws_cf_axi, visDefault.caxis_img), title('SWS axi'), colormap("turbo"), colorbar
subplot(2,2,3), imagesc(sws_cf_both, visDefault.caxis_img), title('SWS both'), colormap("turbo"), colorbar
subplot(2,2,4), imagesc(sws_cf_1d, visDefault.caxis_img), title('SWS 1d'), colormap("turbo"), colorbar


%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizCF = Visualizer_v2(visDefault);

vizCF = vizCF.setROI(xdim, ydim, bigImg(sws_cf_both, Bmode));

titleName = strcat(methodName, ' ID-', idName);
vizCF = vizCF.setTitle(titleName);
vizCF = vizCF.setUnits('mm');

vizCF.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%

%% WAVE APPROXIMATION (MAOW)
if methodMAOW
methodName = 'MAOW';
window = 85; 
w_kernel = [window, window];

% pv_field = Wave_z;  
pv_field = angle(exp_phase);
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

sws_maow = sws_generator(mirror_frame,w_kernel,freq,2,dinf,og_size,10,5);
sws_maow_big = bigImg(sws_maow, Bmode); % resize to Bmode size

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizMAOW = Visualizer_v2(visDefault);

vizMAOW = vizMAOW.setROI(xdim, ydim, sws_maow_big);

titleName = strcat(methodName, ' ID-', idName);
vizMAOW = vizMAOW.setTitle(titleName);
vizMAOW = vizMAOW.setUnits('mm');

vizMAOW.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
end
%% PHASE GRADIENT LEAST SQUARES ORIGINAL PAPER
if methodPG_LS
methodName = 'PG-LS';
stride = 2; % Stride for window-based method
window = 21; 
w_kernel = [window, window];

% pv_field = Wave_z;  
pv_field = exp_phase;
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

constant = 0.33;
[grad_z,grad_x,k,sws_pg_ls] = phase_estimator_lsq(mirror_frame, w_kernel,freq,dinf,og_size,constant);

sws_pg_ls_big = bigImg(sws_pg_ls, Bmode); % resize to Bmode size 

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizPGLS = Visualizer_v2(visDefault);

vizPGLS = vizPGLS.setROI(xdim, ydim, sws_pg_ls_big);

titleName = strcat(methodName, ' ID-', idName);
vizPGLS = vizPGLS.setTitle(titleName);
vizPGLS = vizPGLS.setUnits('mm');

vizPGLS.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
end
%% PHASE GRADIENT WITH TOTAL VARIATION (PG-TV)
if methodPG_TV
methodName = 'PG-TV';

stride = 2; % Stride for window-based method
window = 41; 
w_kernel = [window, window];

% pv_field = Wave_z;  
pv_field = exp_phase;
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

% L2 NORM PG  
[grad_l2, size_out] = pg_norm(mirror_frame, w_kernel, dinf, og_size, stride);
sws_pg = (2*pi*freq)./grad_l2;
sws_pg_big = bigImg(sws_pg, Bmode); % resize to Bmode size

% TV Denoising
mu = 10^2;
tol = 1e-4;
M = size_out(1); N = size_out(2);
[b_opt] = IRLS_TV(grad_l2(:),speye(M*N),mu,M,N,tol,ones(size(M*N)),ones(M*N,1));

grad_l2_tv = reshape(b_opt, size(grad_l2));
sws_pg_tv = (2*pi*freq)./grad_l2_tv;

sws_pg_tv_big = bigImg(sws_pg_tv, Bmode); % resize to Bmode size

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizTV = Visualizer_v2(visDefault);

vizTV = vizTV.setROI(xdim, ydim, sws_pg_tv_big);

titleName = strcat(methodName, ' ID-', idName);
vizTV = vizTV.setTitle(titleName);
vizTV = vizTV.setUnits('mm');

vizTV.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
end
%%
%% NEW VISUALIZATION
visDefault.x = xdim;                % Lateral coordinates
visDefault.z = ydim;                % Axial coordinates
visDefault.BmodeFull = Bmode;       % B-mode image data
visDefault.caxis_bmode = [-60 0];   % Color axis limits for B-mode
visDefault.caxis_img = [0 4.5];     % Color axis limits for color (SWS)
visDefault.fact_transparency = 0.55; % Example transparency factor

vizTV = Visualizer_v2(visDefault);
vizTV = vizTV.setROI(xdim, ydim, sws_pg_tv_big);
titleName = strcat(methodName, ' ID-', idName);
vizTV = vizTV.setTitle(titleName);
vizTV = vizTV.setUnits('mm');
vizTV.visualize();