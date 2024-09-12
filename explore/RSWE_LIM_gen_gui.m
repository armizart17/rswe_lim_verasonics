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

% Wait for the user to load the file before proceeding
typName = matFileLoader; % NEW ONE

% Check if a file was successfully loaded
if isempty(typName)
    error('No file was loaded!');
end

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

%% CALCULATION OF WAVELENGTH
sws_phantom = 4.5; % asume 4.5m/s
[wvlength, wvnumber] = calc_wvlength_k(sws_phantom, freq, dinf);


%% VISUALIZATION SWS IMAGES DEFAULT

% Define default Visualization specs
visDefault.x = xdim;                % Lateral coordinates
visDefault.z = ydim;                % Axial coordinates
visDefault.BmodeFull = Bmode;       % B-mode image data
visDefault.caxis_bmode = [-60 0];   % Color axis limits for B-mode
visDefault.caxis_img = [0 4.5];     % Color axis limits for color (SWS)
visDefault.fact_transparency = 0.6; % Example transparency factor

%% CURVE FITTING (CF) "VERSION LIM"
methodName = 'CF-LIM';

% Fix to 1.5 wv 
factCF = 1.5;
win(1) = round(factCF*wvlength.pix_axi); 
win(2) = round(factCF*wvlength.pix_lat); 
% Make it odd by adding 1 if it's even
if mod(win(1), 2) == 0
   win(1) = win(1) + 1; 
end
if mod(win(2), 2) == 0
   win(2) = win(2) + 1; 
end

correc = xcorr(ones(win(1), win(2))); % it is default 
dx = dinf.dx; dz = dinf.dz;

pv_field = Frames0;  
og_size = size(pv_field);
mirror_pv = padarray(pv_field, (win-1)/2, 'symmetric');

tic 
[Kx,Kz,Rx,Rz] = sws_estimation_curve_fitting(mirror_pv, win, dx , dz, correc);
tt = toc;
fprintf('Time passed CF %.4f\n', tt)

K_tot = 0.5*(Kx + Kz);
sws_cf = real(2*pi*freq/K_tot);
sws_cf_big = bigImg(sws_cf, Bmode); % resize to Bmode size

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizCF = Visualizer(visDefault);

vizCF = vizCF.setROI(xdim, ydim, sws_cf_big);

titleName = strcat(methodName, ' ID-', idName);
vizCF = vizCF.setTitle(titleName);
vizCF = vizCF.setUnits('mm');

vizCF.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%

