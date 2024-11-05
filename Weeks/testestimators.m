win = [101 101];
correc = xcorr2(ones(win(1), win(2))); % it is default 
dx = dinf.dx; dz = dinf.dz;

pv_field = Frames0;

pv_field = u_Filtered;
og_size = size(pv_field);
mirror_pv = padarray(pv_field, (win-1)/2, 'symmetric','both');

tic 
% [Kx,Kz,Rx,Rz] = sws_estimation_curve_fitting(mirror_pv, win, dx , dz, correc);
% [Kx,Kz,Rx,Rz] = sws_estimation_cf(mirror_pv, win, dx, dz, correc);

% [Kx,Kz,Rx,Rz, K1d, R1d] = sws_estimation_cf(mirror_pv, win, dx, dz, correc); % GILMER
[Kx,Kz,Rx,Rz,K1d,R1d] = sws_estimation_cf_fast(mirror_pv, win, dx, dz, correc, og_size); % EMZ

tt = toc;
fprintf('Time passed CF %.4f\n', tt)

figure, 
sgtitle(['Win = ', num2str(win)])
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
subplot(2,2,1), imagesc(sws_cf_lat, visDefault.caxis_img), title('SWS lat'), colormap("jet"), colorbar
subplot(2,2,2), imagesc(sws_cf_axi, visDefault.caxis_img), title('SWS axi'), colormap("jet"), colorbar
subplot(2,2,3), imagesc(sws_cf_both, visDefault.caxis_img), title('SWS both'), colormap("jet"), colorbar
subplot(2,2,4), imagesc(sws_cf_1d, visDefault.caxis_img), title('SWS 1d'), colormap("jet"), colorbar

%%
%% WAVE APPROXIMATION (MAOW)
methodName = 'MAOW';
win = [71 71];

pv_field = Frames0;  
og_size = size(pv_field);
mirror_pv = padarray(pv_field, (win-1)/2, 'symmetric','both');

[sws_maow,kz_maow, kx_maow] = sws_generator(mirror_pv,win,freq,2,dinf,og_size,10,5);

figure, 
subplot(1,2,1), imagesc(kz_maow), title('Kz MAOW'), colorbar
subplot(1,2,2), imagesc(kx_maow), title('Kx MAOW'), colorbar


sws_maow_big = bigImg(sws_maow, Bmode); % resize to Bmode size

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizMAOW = Visualizer(visDefault);

vizMAOW = vizMAOW.setROI(xdim, ydim, sws_maow_big);

titleName = strcat(methodName, ' ID-', idName);
vizMAOW = vizMAOW.setTitle(titleName);
vizMAOW = vizMAOW.setUnits('mm');

vizMAOW.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%

%%
% TEST PG-NORM
methodName = 'PG-L2';
stride = 2; % Stride for window-based method
% window = 31;

% Fix to 0.75 wv 
factPG = 0.1;
window = round(factPG*wvlength.pix_axi); 
% Make it odd by adding 1 if it's even
if mod(window, 2) == 0
   window = window + 1; 
end
w_kernel = [window, window];

pv_field = Frames0;  

pv_field = u_Filtered; % &FGD EMZ
 
 
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


pv_field = u_Filtered; % &FGD EMZ
 


og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

constant = 0.1;
tic;
[grad_z,grad_x,k,sws_pg_ls] = phase_estimator_lsq(mirror_frame, w_kernel,freq,dinf,og_size,constant);
tt = toc;
fprintf('Time passed for %s: %.4f\n', methodName, tt);


sws_pg_z = 2*pi*freq./grad_z;
sws_pg_x = 2*pi*freq./grad_x;

figure, 
subplot(1,2,1), imagesc(sws_pg_z, [0 5]), title('Axi'), colorbar, colormap("jet");
subplot(1,2,2), imagesc(sws_pg_x, [0 5]), title('Lat'), colorbar, colormap("jet");

sws_pg_ls_big = bigImg(sws_pg_ls, Bmode); % resize to Bmode size 

%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%
vizPGLS = Visualizer(visDefault);

vizPGLS = vizPGLS.setROI(xdim, ydim, sws_pg_ls_big);

titleName = strcat(methodName, ' ID-', idName);
vizPGLS = vizPGLS.setTitle(titleName);
vizPGLS = vizPGLS.setUnits('mm');

vizPGLS.visualize(); % This will plot 
%%%%%%%%%%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%%%%%%%%%%%%