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

addpath(genpath(pwd));
%% Setup Directories and FolderName
clear all, clc
close all;

baseDir = 'D:\emirandaz\rswe\data_lim';     %%%% *CHANGE IT IF NEED IT %%%%
folderAcqName = '05-09' ;                   %%%% *CHANGE IT IF NEED IT %%%%

dataDir = fullfile(baseDir, folderAcqName);

fprintf('==============================================================\n');
fprintf('Data Directory Path: %s\n', dataDir);
fprintf('==============================================================\n');

% outDir = ; % for Results 
%% Setup TypicalName of Data

typName = input('Enter full name of data (i.e. RSWE_2lay_500Hz): ', 's');

% typName = "RSWE_"+idName+".mat"; % old version

pos_ = strfind(typName, '_');
posHz = strfind(typName, 'Hz');

idName = typName(pos_(1)+1 : posHz+1);

fprintf('==============================================================\n');
fprintf('Data Full Name: %s\n', typName)
fprintf('Data ID Name: %s\n', idName)
fprintf('==============================================================\n');

% Extract the four digits before "Hz"
freq = str2double( typName(posHz-4:posHz-1) ); % frequency usually represented by 4digits 
%% LOAD DATA AND BMODE VALIDATION
load(fullfile(dataDir, typName +".mat"));

dinf.fc = Trans.frequency*1e6;
dinf.c0 = 1540; dinf.wl = dinf.c0/dinf.fc;
dinf.dz = PData(1).PDelta(3)*dinf.wl;
dinf.dx = PData(1).PDelta(1)*dinf.wl;
dinf.samplesRF = 2048; dinf.PRFe = 10*500;
dinf.Tprf = 1/dinf.PRFe; dinf.Troundtrip = (dinf.samplesRF/2)*(1/dinf.fc);
dinf.maxDepth = dinf.Troundtrip*dinf.c0/2;
dinf.fs = 5000; fs = 5000;dinf.offset_y = 0; 

IQFrame = IQData1(:,:,end); % last frame
Bmode = 20*log10(abs(IQFrame(:,:,1)));
Bmode = Bmode - max(Bmode(:));

dinf.num_angles = 1;
dinf.offset_z = dinf.dz*size(IQFrame,1)/2;

xdim = linspace(-dinf.dx*size(Bmode,2)/2,dinf.dx*size(Bmode,2)/2,size(Bmode,2)); 
ydim = linspace(0,dinf.dz*size(Bmode,1),size(Bmode,1));

figure, 
imagesc(xdim*1e2, ydim*1e2, Bmode), 
colormap('gray'), 
xlabel('\bfLateral [cm]'),
ylabel('\bfAxial [cm]');
title(['B-mode ID: ', idName], 'Interpreter', 'none', 'FontWeight', 'bold')


%% PARTICLE VELOCITY ESTIMATION 

[u,dinf] = pv_cal(IQData1,dinf, dinf.num_angles); 
u = signal_period(freq, dinf.PRFe, u);   

% Temperal filtering process, a bandpass FIR filter is used around 
% +- 20 Hz the vibration frequency
cs_min = 0.7; % [m/s]
cs_max = 5;   % [m/s]
f_tol = 50;   % [Hz] tolerance for +/-2*f_tol
% This function uses spatial_fil_phase_extrac inside


[u_new, Frames0, Frames1] = u_filt(u, freq, f_tol, dinf, cs_min, cs_max);  
% u_new: video
% Frames0: complex PV
% Frames1: onlyPhase PV


% CHECK PLOT FRAMES0 and FRAMES 1
% figure(), 
% sgtitle(['Frames0 ID: ', idName], 'Interpreter', 'none')
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
% subplot(121), 
% imagesc(xdim*1e2, ydim*1e2,abs(Frames0)), title('Magn'), axis("tight");
% xlabel('Lateral [cm]'), ylabel('Axial [cm]')
% % set(gca, "clim", [1, 5])
% subplot(122), 
% imagesc(xdim*1e2, ydim*1e2, angle(Frames0)), title('Phase'), axis("tight"), colormap('default')
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
methodName = 'CF';
window = 21; 
w_kernel = [window, window];

pv_field = Frames0;  
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

[k_z,R_ax,k_x,R_lat,k,sws_cf] = theoretical_fitting(mirror_frame,w_kernel,freq,dinf,og_size);
sws_cf_big = bigImg(sws_cf, Bmode); % resize to Bmode size

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizCF = Visualizer(visDefault);

vizCF = vizCF.setROI(xdim, ydim, sws_cf_big);

titleName = strcat(methodName, ' ID-', idName);
vizCF = vizCF.setTitle(titleName);
vizCF = vizCF.setUnits('mm');

vizCF.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%

%% WAVE APPROXIMATION (MAOW)
methodName = 'MAOW';
window = 25; 
w_kernel = [window, window];

pv_field = Frames0;  
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

sws_maow = sws_generator(mirror_frame,w_kernel,freq,2,dinf,og_size,10,5);
sws_maow_big = bigImg(sws_pg_tv, Bmode); % resize to Bmode size

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizMAOW = Visualizer(visDefault);

vizMAOW = vizMAOW.setROI(xdim, ydim, sws_maow_big);

titleName = strcat(methodName, ' ID-', idName);
vizMAOW = vizMAOW.setTitle(titleName);
vizMAOW = vizMAOW.setUnits('mm');

vizMAOW.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%

%% PHASE GRADIENT LEAST SQUARES ORIGINAL PAPER
methodName = 'PG-LS';
stride = 2; % Stride for window-based method
window = 15; 
w_kernel = [window, window];

pv_field = Frames0;  
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

constant = 0.33;
[grad_z,grad_x,k,sws_pg_ls] = phase_estimator_lsq(mirror_frame, w_kernel,freq,dinf,og_size,constant);

sws_pg_ls_big = bigImg(sws_pg_ls, Bmode); % resize to Bmode size 

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizPGLS = Visualizer(visDefault);

vizPGLS = vizPGLS.setROI(xdim, ydim, sws_pg_ls_big);

titleName = strcat(methodName, ' ID-', idName);
vizPGLS = vizPGLS.setTitle(titleName);
vizPGLS = vizPGLS.setUnits('mm');

vizPGLS.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%

%% PHASE GRADIENT WITH TOTAL VARIATION (PG-TV)
methodName = 'PG-TV';

stride = 2; % Stride for window-based method
window = 41; 
w_kernel = [window, window];

pv_field = Frames0;  
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

% L2 NORM PG  
[grad_l2, size_out] = pg_norm(mirror_frame, w_kernel, dinf, og_size, stride);
sws_pg = (2*pi*freq)./grad_l2;
sws_pg_big = bigImg(sws_pg, Bmode); % resize to Bmode size

% TV Denoising
mu = 10^4;
tol = 1e-4;
M = size_out(1); N = size_out(2);
[b_opt] = IRLS_TV_simple(grad_l2(:),speye(M*N),mu,M,N,tol,ones(size(M*N)),ones(M*N,1));

grad_l2_tv = reshape(b_opt, size(grad_l2));
sws_pg_tv = (2*pi*freq)./grad_l2_tv;

sws_pg_tv_big = bigImg(sws_pg_tv, Bmode); % resize to Bmode size

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizTV = Visualizer(visDefault);

vizTV = vizTV.setROI(xdim, ydim, sws_pg_tv_big);

titleName = strcat(methodName, ' ID-', idName);
vizTV = vizTV.setTitle(titleName);
vizTV = vizTV.setUnits('mm');

vizTV.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%

%%
