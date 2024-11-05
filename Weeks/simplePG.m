% SIMPLE SCRIPT 4 GFB BY EMZ
% Sacar functiones de las otras carpetas
% phase_estimator_lsq ...... /SWS_PhaseGradient
% pg_norm ...... /SWS_PhaseGradient
% IRLS_TV_simple .... /SWS_PhaseGradient

% Cosas cambiar
% window, Frames0
% stride = 2


%% PHASE GRADIENT LEAST SQUARES ORIGINAL PAPER
methodName = 'PG-LS';
stride = 2; % Stride for window-based method


% window = 31;
w_kernel = [window, window];

pv_field = Frames0;  
og_size = size(pv_field);
mirror_frame = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

constant = 0.1;
tic;
[grad_z,grad_x,k,sws_pg_ls] = phase_estimator_lsq(mirror_frame, w_kernel,freq,dinf,og_size,constant);
tt = toc;
fprintf('Time passed for %s: %.4f\n', methodName, tt);

figure, 
imagesc(sws_pg_ls)
title(methodName)

%%
%%
% REAL PHANTOM
% pv_field = Frames0; 

x = xdim;
z = ydim;

cs_min = 1; cs_max = 8;
calc_wvlength_k(cs_min, freq, dinf);
calc_wvlength_k(cs_max, freq, dinf);
% [u_filt_v2] = fun_JO_v2(u, freq, dinf, cs_min, cs_max);

% [M_filtered, u_filt_v2] = fun_JO_v1(u, freq, dinf);
[M_filtered, u_filt_v2] = fun_JO_v2(u, freq, dinf, cs_min, cs_max);

pv_field = u_filt_v2;

% SIMULATION
% pv_field = pv_complexZ;
% dinf.dx = min(diff(x));
% dinf.dz = min(diff(z));

[Ny, Nx, ~] = size(pv_field);

% Define the spatial sampling intervals (assuming equal spacing)
dx = dinf.dx;  % Example spatial resolution in meters
dy = dinf.dz;  % Assume square grid for simplicity

% Compute the wavenumber ranges
kx = (-Nx/2:Nx/2-1) * (2*pi / (Nx * dx));
ky = (-Ny/2:Ny/2-1) * (2*pi / (Ny * dy));

% Create the wavenumber grid
[KX, KY] = meshgrid(kx, ky); 


% Perform the 2D FFT and shift the zero-frequency component to the center
fft_velocity = fft2(pv_field);
fft_shifted = fftshift(fft_velocity);
V_k = abs(fft_shifted);


figure, 
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.8 0.55])
sgtitle(['Frames1 ID: ', idName], 'Interpreter', 'none')

subplot(1,2,2)
imagesc( kx, ky, log(1+V_k)  ), grid on;
colorbar, 
title('K map')
xlabel('Kx (rad/m)'), ylabel('Ky (rad/m)')
xlim([-8000 8000])
ylim([-8000 8000])

subplot(1,2,1)
imagesc( x*1e3, z*1e3, angle(pv_field)  );
colorbar, 
title(['k = ', num2str(wvnumber)])
xlabel('Lateral (mm)'), ylabel('Depth (mm)')
% xlim([-5000 5000])
% ylim([-5000 5000])

%% PHASE GRADIENT L2-norm
methodName = 'PG-L2';
stride = 2; % Stride for window-based method
window = 15;

w_kernel = [window, window];


% [Frames0] = fun_JO_v2(u, freq, dinf);

% pv_field = Frames0;
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

sws_big = bigImg(sws_pg, Bmode);

mean_sws = mean(sws_big(:));
std_sws = std(sws_big(:));

pmSymbol = char(177); %+/-
title_str = sprintf('%s: %.2f %c %.2f', methodName, mean_sws, pmSymbol, std_sws);
            

visDefault.caxis_bmode = [-60 0 ];
visDefault.caxis_img = [0 4];

figure,
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.8 0.55])
sgtitle('\bf L2-Phase Gradient ||\nabla\phi||_2')

subplot(1,2,1)
imagesc(xdim*1e3, ydim*1e3, sws_big, visDefault.caxis_img)
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title(title_str);
axis("tight")
colorbar;

fact_transparency = 0.5;
subplot(1,2,2)
[hF,hB,hColor] = imOverlayInterp(Bmode,sws_big,visDefault.caxis_bmode, visDefault.caxis_img, ...
                fact_transparency,xdim*1e3,ydim*1e3,sws_big,xdim*1e3,ydim*1e3);
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title(title_str);

% PHASE GRADIENT WITH TOTAL VARIATION (PG-TV)
methodName = 'PG-TV';

% Fix to 0.5 wv 
w_kernel = [window, window];

% pv_field = Frames0;  
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


sws_big = bigImg(sws_pg_tv, Bmode);
mean_sws = mean(sws_big(:));
std_sws = std(sws_big(:));

pmSymbol = char(177); %+/-
title_str = sprintf('%s: %.2f %c %.2f', methodName, mean_sws, pmSymbol, std_sws);
            

figure,
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.8 0.55])
sgtitle('\bf L2-PG-TV ||\nabla\phi||_2')

subplot(1,2,1)
imagesc(xdim*1e3, ydim*1e3, sws_big, visDefault.caxis_img)
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title(title_str);
axis("tight")
colorbar;

fact_transparency = 0.5;
subplot(1,2,2)
[hF,hB,hColor] = imOverlayInterp(Bmode,sws_big,visDefault.caxis_bmode, visDefault.caxis_img, ...
                fact_transparency,xdim*1e3,ydim*1e3,sws_big,xdim*1e3,ydim*1e3);
xlabel('Lateral [mm]'), ylabel('Depth [mm]')
title(title_str);




%%
figure, 
plot(ff, cs_min_ff, 'DisplayName','cmin')
hold on;
plot(ff, cs_max_ff, 'DisplayName','cmax')
hold off;