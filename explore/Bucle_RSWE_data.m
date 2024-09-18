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

%%
%% Setup Directories and FolderName
clear all, clc
close all;

baseDir = 'D:\emirandaz\rswe\data_lim\12-09\300_360_400\';     %%%% *CHANGE IT IF NEED IT %%%%
folderAcqName = 'both' ;                   %%%% *CHANGE IT IF NEED IT %%%%

dataDir = fullfile(baseDir, folderAcqName);

figDir = fullfile(dataDir,'fig');
if ~exist("figDir","dir"); mkdir(figDir); end

fprintf('==============================================================\n');
fprintf('Data Directory Path: %s\n', dataDir);
fprintf('==============================================================\n');

% outDir = ; % for Results 
%% Setup TypicalName of Data (this time process all)

data_files = dir([dataDir,'/*.mat']);
num_files = length(data_files);

% num_files = [1];
for ii = 1:num_files

%%% BUCLE %%%
fileName  = data_files(ii).name;  % "typName.mat"
typName = fileName(1:end-4);      % 
%%% BUCLE %%%

% typName = input('Enter full name of data (i.e. RSWE_2lay_500Hz): ', 's');

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
figure(), 
sgtitle(['Frames0 ID: ', idName], 'Interpreter', 'none')
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
subplot(121), 
imagesc(xdim*1e2, ydim*1e2,abs(Frames0)), title('Magn'), axis("tight");
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
% set(gca, "clim", [1, 5])
subplot(122), 
imagesc(xdim*1e2, ydim*1e2, angle(Frames0)), title('Phase'), axis("tight"), colormap('default')
xlabel('Lateral [cm]'), ylabel('Axial [cm]')

nameFig = "pv_" + idName;
saveas(gcf,fullfile(figDir, nameFig + ".png"));

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

visDefault.fact_transparency = 0.6; % Example transparency factor

%% CALCULATION OF WAVELENGTH
% sws_phantom = 5; % 12%
tokens = regexp(idName, '(\d+conc)_\d+Hz', 'tokens');
concentrationPhantom = tokens{1}{1};

switch concentrationPhantom
        case '16conc'
            sws_phantom = 6.2;
            visDefault.caxis_img = [0 8];     % Color axis limits for color (SWS)
        case '12conc'
            sws_phantom = 5;
            visDefault.caxis_img = [0 6.5];     % Color axis limits for color (SWS)
        case '10conc'
            sws_phantom = 4;  % You can also choose to average 4 and 5 here if needed
        case '4conc'
            sws_phantom = 3;
        otherwise
            error('Invalid phantom percentage value.');
end
[wvlength, wvnumber] = calc_wvlength_k(sws_phantom, freq, dinf);

disp(concentrationPhantom);
fprintf('SWS %.2f m/s | Freq %d Hz | Kernel size Ax: %d | La: %d \n', ...
     sws_phantom, freq, round(wvlength.pix_axi), round(wvlength.pix_lat)  );

%% PHASE GRADIENT LEAST SQUARES ORIGINAL PAPER
methodName = 'PG-LS';
stride = 2; % Stride for window-based method
% window = 31;

% Fix to 0.75 wv 
factPG = 0.5;
window = round(factPG*wvlength.pix_axi); 
% Make it odd by adding 1 if it's even
if mod(window, 2) == 0
   window = window + 1; 
end
w_kernel = [window, window];

pv_field = Frames0;  
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

constant = 0.1;
tic;
[grad_z,grad_x,k,sws_pg_ls] = phase_estimator_lsq(mirror_frame, w_kernel,freq,dinf,og_size,constant);
tt = toc;
fprintf('Time passed for %s: %.4f\n', methodName, tt);

sws_pg_ls_big = bigImg(sws_pg_ls, Bmode); % resize to Bmode size 

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizPGLS = Visualizer(visDefault);

vizPGLS = vizPGLS.setROI(xdim, ydim, sws_pg_ls_big);

titleName = strcat(methodName, ' ID-', idName);
vizPGLS = vizPGLS.setTitle(titleName);
vizPGLS = vizPGLS.setUnits('mm');

vizPGLS.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%

nameFig = strcat(methodName, '_',idName);
saveas(gcf,fullfile(figDir, nameFig + ".png"));
save(fullfile(figDir, nameFig)+".mat", "xdim", "ydim", "sws_pg_ls_big", "Bmode");
%% PHASE GRADIENT L2-noem
methodName = 'PG-L2';
stride = 2; % Stride for window-based method
% window = 31;

% Fix to 0.75 wv 
factPG = 0.5;
window = round(factPG*wvlength.pix_axi); 
% Make it odd by adding 1 if it's even
if mod(window, 2) == 0
   window = window + 1; 
end
w_kernel = [window, window];

pv_field = Frames0;  
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

% L2 NORM PG
tic;
[grad_l2, size_out] = pg_norm(mirror_frame, w_kernel, dinf, og_size, stride);
sws_pg = (2*pi*freq)./grad_l2;
tt = toc;
fprintf('Time passed for %s: %.4f\n', methodName, tt);


% POST-PROCESSING 
kernel_post = 7;

% AVERAGE FILTER
% avg_kernel = ones(kernel_post, kernel_post) / kernel_post^2;
% sws_pg = filter2(avg_kernel, sws_pg, 'same');

% MEDIAN FILTER
sws_pg = medfilt2(sws_pg,[kernel_post kernel_post],'symmetric');

sws_pg_big = bigImg(sws_pg, Bmode); % resize to Bmode size

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizPGL2 = Visualizer(visDefault);

vizPGL2 = vizPGL2.setROI(xdim, ydim, sws_pg_big);

titleName = strcat(methodName, ' ID-', idName);
vizPGL2 = vizPGL2.setTitle(titleName);
vizPGL2 = vizPGL2.setUnits('mm');

vizPGL2.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%

nameFig = strcat(methodName, '_',idName);
saveas(gcf,fullfile(figDir, nameFig + ".png"));
save(fullfile(figDir, nameFig)+".mat", "xdim", "ydim", "sws_pg_big", "Bmode");

%% PHASE GRADIENT WITH TOTAL VARIATION (PG-TV)
methodName = 'PG-TV';

% Fix to 0.5 wv 

window = round(factPG*wvlength.pix_axi); 
% Make it odd by adding 1 if it's even
if mod(window, 2) == 0
   window = window + 1; 
end
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

tic;
[b_opt] = IRLS_TV_simple(grad_l2(:),speye(M*N),mu,M,N,tol,ones(size(M*N)),ones(M*N,1));
tt = toc;
fprintf('Time passed for %s: %.4f\n', methodName, tt);

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

nameFig = strcat(methodName, '_',idName);
saveas(gcf,fullfile(figDir, nameFig + ".png"));
save(fullfile(figDir, nameFig)+".mat", "xdim", "ydim", "sws_pg_tv_big", "Bmode");


%% CURVE FITTING (CF) "VERSION LIM"
methodName = 'CF-LIM';

% Fix to 1.1 wv 
factCF = 1.1;
win(1) = round(factCF*wvlength.pix_axi); 
win(2) = round(factCF*wvlength.pix_lat); 
% Make it odd by adding 1 if it's even
if mod(win(1), 2) == 0
   win(1) = win(1) + 1; 
end
if mod(win(2), 2) == 0
   win(2) = win(2) + 1; 
end

correc = xcorr2(ones(win(1), win(2))); % it is default 
dx = dinf.dx; dz = dinf.dz;

pv_field = Frames0;  
og_size = size(pv_field);
mirror_pv = padarray(pv_field, (win-1)/2, 'symmetric','both');

tic 
% [Kx,Kz,Rx,Rz] = sws_estimation_curve_fitting(mirror_pv, win, dx , dz, correc);
[Kx,Kz,Rx,Rz] = sws_estimation_cf(mirror_pv, win, dx, dz, correc);
tt = toc;
fprintf('Time passed for %s: %.4f\n', methodName, tt);

K_tot = 0.5*(Kx + Kz);
sws_cf = real(2*pi*freq./K_tot);
sws_cf_big = bigImg(sws_cf, Bmode); % resize to Bmode size

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizCF = Visualizer(visDefault);

vizCF = vizCF.setROI(xdim, ydim, sws_cf_big);

titleName = strcat(methodName, ' ID-', idName);
vizCF = vizCF.setTitle(titleName);
vizCF = vizCF.setUnits('mm');

vizCF.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%

nameFig = strcat(methodName, '_',idName);
saveas(gcf,fullfile(figDir, nameFig + ".png"));
save(fullfile(figDir, nameFig)+".mat", "xdim", "ydim", "sws_cf_big", "Bmode");

%% CURVE FITTING (CF) "VERSION EMZ"
% methodName = 'CF';
% window = 21; 
% 
% Fix to 1.5 wv 
% factCF = 1.5;
% window = round(factCF*wvlength.pix_axi); 
% Make it odd by adding 1 if it's even
% if mod(window, 2) == 0
%    window = window + 1; 
% end
% 
% w_kernel = [window, window];
% 
% pv_field = Frames0;  
% og_size = size(pv_field);
% mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');
% 
% tic
% [k_z,R_ax,k_x,R_lat,k,sws_cf] = theoretical_fitting(mirror_frame,w_kernel,freq,dinf,og_size);
% sws_cf_big = bigImg(sws_cf, Bmode); % resize to Bmode size
% tt = toc;
% fprintf('Time passed CF %.4f\n', tt)
% %%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
% vizCF = Visualizer(visDefault);
% 
% vizCF = vizCF.setROI(xdim, ydim, sws_cf_big);
% 
% titleName = strcat(methodName, ' ID-', idName);
% vizCF = vizCF.setTitle(titleName);
% vizCF = vizCF.setUnits('mm');
% 
% vizCF.visualize(); % This will plot 
% %%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
% % WAVE APPROXIMATION (MAOW)
% 
% methodName = 'MAOW';
% 
% Fix to 1.5 wv 
% factMaow = 1.5;
% window = round(factMaow*wvlength.pix_axi); 
% Make it odd by adding 1 if it's even
% if mod(window, 2) == 0
%    window = window + 1; 
% end
%  
% w_kernel = [window, window];
% 
% pv_field = Frames0;  
% og_size = size(pv_field);
% mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');
% 
% tic
% sws_maow = sws_generator(mirror_frame,w_kernel,freq,2,dinf,og_size,10,5);
% sws_maow_big = bigImg(sws_maow, Bmode); % resize to Bmode size
% tt = toc;
% fprintf('Time passed MAOW %.4f\n', tt)
% 
% %%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
% vizMAOW = Visualizer(visDefault);
% 
% vizMAOW = vizMAOW.setROI(xdim, ydim, sws_maow_big);
% 
% titleName = strcat(methodName, ' ID-', idName);
% vizMAOW = vizMAOW.setTitle(titleName);
% 
% vizMAOW = vizMAOW.setUnits('mm');
% 
% vizMAOW.visualize(); % This will plot 
% %%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%


close all;
end