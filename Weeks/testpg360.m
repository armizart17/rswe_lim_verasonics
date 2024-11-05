%%


[JO_filtered] = fun_JO_v1(u, freq, dinf);


u_T = permute(u, [2, 1, 3]);
[JO_filtered_T] = fun_JO_v1(u_T, freq, dinf);


figure(), 
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
sgtitle(['JO ID u'], 'Interpreter', 'none')
% sgtitle(['Frames1 ID: ', idName], 'Interpreter', 'none')

subplot(121), 
imagesc(abs(JO_filtered)), title('Mang'), axis("tight");
% set(gca, "clim", [1, 5])
subplot(122), 
imagesc(angle(JO_filtered)), title('Phase'), axis("tight"), colormap('default')

figure(), 
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
sgtitle(['JO^T ID u'], 'Interpreter', 'none')
% sgtitle(['Frames1 ID: ', idName], 'Interpreter', 'none')

subplot(121), 
imagesc(abs(JO_filtered_T)), title('Mang'), axis("tight");
% set(gca, "clim", [1, 5])
subplot(122), 
imagesc(angle(JO_filtered_T)), title('Phase'), axis("tight"), colormap('default')


%%
[uLIM2,dinf] = pv_cal(IQData1,dinf, dinf.num_angles); 

uLIM2 = signal_period(freq, dinf.PRFe, uLIM2);   
cs_min = 0.7; % [m/s]
cs_max = 5;   % [m/s]
f_tol = 50;   % [Hz] tolerance for +/-2*f_tol
[u_new, Frames0, Frames1] = u_filt(uLIM, freq, f_tol, dinf, cs_min, cs_max);  


figure(), 
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
sgtitle(['Frames0 ID'], 'Interpreter', 'none')
% sgtitle(['Frames1 ID: ', idName], 'Interpreter', 'none')

subplot(121), 
imagesc(abs(Frames0)), title('Mang'), axis("tight");
% set(gca, "clim", [1, 5])
subplot(122), 
imagesc(angle(Frames0)), title('Phase'), axis("tight"), colormap('default')


%% TEST DATA JUVENAL 
load('P:\rswe\Data_3_PUCP\CIRS_phantom\L7-4\data_13\data_13.mat')
dinf.offset_y = 0;
freq = freq(1);


% % [uLIM,dinf] = pv_cal(IQData1,dinf, dinf.num_angles); 
% uLIM = u;
% uLIM = signal_period(freq, dinf.PRFe, uLIM);   
% 
% % Temperal filtering process, a bandpass FIR filter is used around 
% % +- 20 Hz the vibration frequency
% cs_min = 0.7; % [m/s]
% cs_max = 5;   % [m/s]
% f_tol = 50;   % [Hz] tolerance for +/-2*f_tol
% % This function uses spatial_fil_phase_extrac inside
% 
% [u_new, Frames0, Frames1] = u_filt(uLIM, freq, f_tol, dinf, cs_min, cs_max);  
% 
% figure(), 
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
% sgtitle(['Frames1 ID'], 'Interpreter', 'none')
% % sgtitle(['Frames1 ID: ', idName], 'Interpreter', 'none')
% 
% subplot(121), 
% imagesc(abs(Frames0)), title('Mang'), axis("tight");
% % set(gca, "clim", [1, 5])
% subplot(122), 
% imagesc(angle(Frames0)), title('Phase'), axis("tight"), colormap('default')


[JO_filtered] = fun_JO_v1(u, freq, dinf);

figure(), 
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
sgtitle(['JO ID u'], 'Interpreter', 'none')
% sgtitle(['Frames1 ID: ', idName], 'Interpreter', 'none')

subplot(121), 
imagesc(abs(JO_filtered)), title('Mang'), axis("tight");
% set(gca, "clim", [1, 5])
subplot(122), 
imagesc(angle(JO_filtered)), title('Phase'), axis("tight"), colormap('default')
% 
% [JO_filtered] = fun_JO_v1(uLIM, freq, dinf);
% 
% figure(), 
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
% sgtitle(['JO ID uLIM'], 'Interpreter', 'none')
% % sgtitle(['Frames1 ID: ', idName], 'Interpreter', 'none')
% 
% subplot(121), 
% imagesc(abs(JO_filtered)), title('Mang'), axis("tight");
% % set(gca, "clim", [1, 5])
% subplot(122), 
% imagesc(angle(JO_filtered)), title('Phase'), axis("tight"), colormap('default')


Bmode = 20*log10(abs(IQBmodeData(:,:,1)));
Bmode = Bmode - max(Bmode(:));

xdim = linspace(-dinf.dx*size(Bmode,2)/2,dinf.dx*size(Bmode,2)/2,size(Bmode,2)); 
ydim = linspace(0,dinf.dz*size(Bmode,1),size(Bmode,1));

%% CURVE FITTING
win = [31 15];
correc = xcorr2(ones(win(1), win(2))); % it is default 
dx = dinf.dx; dz = dinf.dz;

pv_field = angle(JO_filtered);

og_size = size(pv_field);
mirror_pv = padarray(pv_field, (win-1)/2, 'symmetric','both');

tic 
% [Kx,Kz,Rx,Rz] = sws_estimation_curve_fitting(mirror_pv, win, dx , dz, correc);
% [Kx,Kz,Rx,Rz] = sws_estimation_cf(mirror_pv, win, dx, dz, correc);

% [Kx,Kz,Rx,Rz, K1d, R1d] = sws_estimation_cf(mirror_pv, win, dx, dz, correc); % GILMER
[Kx,Kz,Rx,Rz, K1d, R1d] = sws_estimation_cf_fast(mirror_pv, win, dx, dz, correc, og_size); % EMZ
wvnumber = 1000;
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
subplot(2,2,1), imagesc(sws_cf_lat, [0 4]), title('SWS lat'), colormap("jet"), colorbar
subplot(2,2,2), imagesc(sws_cf_axi, [0 4]), title('SWS axi'), colormap("jet"), colorbar
subplot(2,2,3), imagesc(sws_cf_both, [0 4]), title('SWS both'), colormap("jet"), colorbar
subplot(2,2,4), imagesc(sws_cf_1d, [0 4]), title('SWS 1d'), colormap("jet"), colorbar

%% VISUALIZATION SIMPLE
visDefault.caxis_bmode = [-60 0 ];
visDefault.caxis_img = [0 4];
% VISUALIZATION IMOVERLAY

sws_cf_big = bigImg(sws_cf_1d, Bmode);
figure, 
[hF,hB,hColor] = imOverlayInterp(Bmode,sws_cf_big,visDefault.caxis_bmode, visDefault.caxis_img, ...
                0.2,xdim*1e3,ydim*1e3,sws_cf_big,xdim*1e3,ydim*1e3);
title('CF')


%% VISUALIZATION SWS IMAGES DEFAULT

win = [101 101];
og_size = size(pv_complexZ);
mirror_frame = padarray(pv_complexZ, (win-1)/2, 'symmetric');
correc = xcorr2(ones(win(1), win(2))); % it is default 

dinf_sim.dx = min(diff(x));
dinf_sim.dz = min(diff(z));

tic
[Kx,Kz,Rx,Rz, K1d, R1d] = sws_estimation_cf_fast(mirror_frame, win, dinf_sim.dx, dinf_sim.dz, correc, og_size); % EMZ
tt = toc;
fprintf('Time passed CF LIM %.4f\n', tt)
sws_cf_lat = 2*pi*freq./Kx;
sws_cf_axi = 2*pi*freq./Kz;
sws_cf_both = 2*pi*freq./( 0.5*(Kz+Kx)  );
sws_cf_1d = 2*pi*freq./K1d;

figure, 
sgtitle(['Win = ', num2str(win)])
subplot(2,2,1), imagesc(x,z, sws_cf_lat, [0 5]), title('SWS lat'), colormap("jet"), colorbar
subplot(2,2,2), imagesc(x,z, sws_cf_axi, [0 5]), title('SWS axi'), colormap("jet"), colorbar
subplot(2,2,3), imagesc(x,z, sws_cf_both, [0 5]), title('SWS both'), colormap("jet"), colorbar
subplot(2,2,4), imagesc(x,z, sws_cf_1d, [0 5]), title('SWS 1d'), colormap("jet"), colorbar

%
load('C:\Users\emirandaz\OneDrive - pucp.pe\Documentos\MATLAB\TESIS_MPSID\R-SWE\phasegradientRSWE\dataold\Data900Hz-10000wvs\R-FIELD_inc_1.mat')

win = [101 101];
og_size = size(pv_complexZ);
mirror_frame = padarray(pv_complexZ, (win-1)/2, 'symmetric');

% FUNCION SWS MAP (Con ajuste de curva)
dinf_sim.dx = min(diff(x));
dinf_sim.dz = min(diff(z));
tic;
[k_z,R_ax,k_x,R_lat,k,sws_matrix] = theoretical_fitting(mirror_frame,win, freq, dinf_sim, og_size);
tt = toc;
fprintf('Time passed CF EMZ %.4f\n', tt)

figure, 
imagesc(x, z, sws_matrix, [0 5]), colormap("jet"), colorbar
title('CF EMZ')

%%

sws_matrix_maow = sws_generator(mirror_frame,win,freq,2,dinf,og_size,10,5);
figure, 
imagesc(sws_matrix_maow), colormap jet
%%
% %% PHASE GRADIENT L2-norm

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

% L2 NORM PG
pv_field = Frames0;  
og_size = size(pv_field);
mirror_pv = padarray(pv_field,[(window-1)/2 (window-1)/2],'symmetric');

tic;
[grad_l2, size_out] = pg_norm(mirror_pv, w_kernel, dinf, og_size, stride);
clear sws_pg
sws_pg = (2*pi*freq)./grad_l2;
tt = toc;
fprintf('Time passed for %s: %.4f\n', methodName, tt);

figure, imagesc(xdim*1e2, ydim*1e2, sws_pg), colormap("jet")

%% SPLIT FUNCTION 

st = stride;
res_z = dinf.dz; % Axial resolution
res_x = dinf.dx; % Lateral resolution

 % Axis of the kernels
z_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(1))*res_z; 
x_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(2))*res_x; 

% Initializing vectors and matrixes    
angle_u = angle(mirror_pv);

% HELP FOR DIMENSIONS %% % THEORY
size_mirror = size(mirror_pv); % ogsize + w_kernel - 1; %  = size (u)
numkernels = floor( (size_mirror - w_kernel)./st + 1 ); % numKernels
size_out = floor( (og_size - 1)./st + 1 );
	     
cont_kernel = 1; 
angle_z = unwrap(angle_u,[],1);
angle_x = unwrap(angle_u,[],2);

figure, 
subplot(1,3,1), imagesc(unwrap(angle_u)), title('angle') ;
subplot(1,3,2), imagesc(angle_z), title('angle z') ;
subplot(1,3,3), imagesc(angle_x), title('angle x');

grad_abs = zeros(size_out);
grad_z = zeros(size_out);
grad_x = zeros(size_out);

% NORMAL LOOP
i_out = 1;
for ii = 1:st:og_size(1)

    j_out = 1;

    for jj = 1:st:og_size(2) %% for faster computing pararell toolbox

        area_z = angle_z(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
        area_x = angle_x(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel

        diff_area_z = diff(area_z, 1, 1) ./ res_z;
        diff_area_x = diff(area_x, 1, 2) ./ res_x;

        diff_area_z = [diff_area_z; diff_area_z(end, :)];  % replicate the last row
        diff_area_x = [diff_area_x, diff_area_x(:, end)];  % replicate the last column

        % figure, 
        % subplot(2,2,1), imagesc(area_z), title('area_z')
        % subplot(2,2,2), imagesc(area_x), title('area_x')
        % subplot(2,2,3), imagesc(diff_area_z)
        % subplot(2,2,4), imagesc(diff_area_x)

        val_abs = sqrt(diff_area_z.^2 + diff_area_x.^2);

        grad_x (i_out, j_out) = mean(diff_area_x, 'all');
        grad_z (i_out, j_out) = mean(diff_area_z, 'all');

        grad_abs(i_out,j_out) = mean(val_abs, 'all');

        j_out = j_out + 1;

    end
   i_out = i_out + 1;
end

figure, 
subplot(1,3,1), imagesc((grad_x)), title('kx'), colorbar, colormap turbo
subplot(1,3,2), imagesc((grad_z)), title('kz'), colorbar, colormap turbo

grad_xz = sqrt( (grad_x).^2 + (grad_z).^2 );
subplot(1,3,3), imagesc(grad_xz), title('l2'), colorbar, colormap turbo

%%
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

%%

constant = 1/10;
tic;
[grad_z,grad_x,k,sws_pg_ls] = phase_estimator_lsq(mirror_pv, w_kernel,freq,dinf,og_size,constant);
tt = toc;
fprintf('Time passed for %s: %.4f\n', methodName, tt);

sws_pg_ls_big = bigImg(sws_pg_ls, Bmode); % resize to Bmode size 

figure, 
subplot(1,3,1), imagesc(abs(grad_x)), title('kx'), colorbar, colormap turbo
subplot(1,3,2), imagesc(abs(grad_z)), title('kz'), colorbar, colormap turbo

grad_xz = sqrt( abs(grad_x).^2 + abs(grad_z).^2 )
subplot(1,3,3), imagesc(grad_xz), title('l2'), colorbar, colormap turbo

figure,  
imagesc(xdim, xdim, sws_pg_ls_big), colormap("jet"), colorbar

%%
%% K estimation using the R-SWE approach

Bmode = 20*log10(abs(IQBmodeData(:,:,1)));
Bmode = Bmode - max(Bmode(:));

xdim = linspace(-dinf.dx*size(Bmode,2)/2,dinf.dx*size(Bmode,2)/2,size(Bmode,2)); 
ydim = linspace(0,dinf.dz*size(Bmode,1),size(Bmode,1));
% this is different for simulated data, to increase resolution

win = [61 61]
inc_res_z = 1;
inc_res_x = 1;
M = win(1); N = win(2);

x_v = xdim;
z_v = ydim;

M_filtered1 = (Frames1);

% Create zero crossing
[a,b,~] = size(M_filtered1);
mask_new0 = zeros([win,0]+[a,b,0]+ [0,0,length(freq)]); mask_new1 = mask_new0;
mask_new0(round(M/2):end-1-floor(M/2),round(N/2):end-1-floor(N/2),:) = (M_filtered1);
[a2,b2,c2] = size(mask_new0);

% Define delta shift in pixels and constants
delta = 3;
c_z = 10;  
c_x = 5;
% c_x = 15/(1*pi);
% c_z = sqrt(5/sqrt(1*pi)); c_x = c_z;

% generate the loop for the autocorrelation
al_pos_z = round(M/2):a2-round(M/2);          % allow locations to use in axial direction
ev_index_z = 1:2:floor(length(al_pos_z)/2)*2; % even index vector
ev_al_pos_z = al_pos_z(ev_index_z);           % even positions from allow locations in axial direction
search_area_z = -round(M/2)+1:round(M/2)-1;   % search index from kernel:axial

al_pos_x = round(N/2):b2-round(N/2);          % allow locations to use in lateral direction
ev_index_x = 1:2:floor(length(al_pos_x)/2)*2; % even index vector
ev_al_pos_x = al_pos_x(ev_index_x);           % even positions from allow locations in lateral direction
search_area_x = -round(N/2)+1:round(N/2)-1;   % search index from kernel:lateral

% Define varibales T
MM = length(search_area_z);
NN = length(search_area_x);
T = xcorr2(ones(MM,NN)); % Triangulo function matrix
T_0 = T(floor((2*MM+1)/2),floor((2*NN+1)/2));
T_delta_x = T(floor((2*MM+1)/2),floor((2*NN+1)/2)+delta);
T_delta_z = T(floor((2*MM+1)/2)+delta,floor((2*NN+1)/2));
% Define the x and z dimension in meters as a function of autocorrelation
% size
z_cm = linspace(-(dinf.dz/inc_res_z)*(2*MM-1)/2,(dinf.dz/inc_res_z)*(2*MM-1)/2,(2*MM-1));
x_cm = linspace(-(dinf.dx/inc_res_x)*(2*NN-1)/2,(dinf.dx/inc_res_x)*(2*NN-1)/2,(2*NN-1));
% Define delta_x and delta_z in meters
delta_z = z_cm(floor((2*MM+1)/2)+delta);
delta_x = x_cm(floor((2*NN+1)/2)+delta);
kz = zeros(length(ev_al_pos_z),length(ev_al_pos_x));
kx = kz;
cs = zeros(length(ev_al_pos_z),length(ev_al_pos_x),length(freq));

    % for kk=1:length(ev_al_pos_z)
    %     s_2D = mask_new0(ev_al_pos_z(kk)+search_area_z,:);
    %     parfor ll=1:length(ev_al_pos_x)
    %         A_c = xcorr2(s_2D(:,ev_al_pos_x(ll)+search_area_x)); % autocorrelation matrix
    %         % Define zp and xp vectors and since z(0)= x(0), then they are defined as
    %         % ZX_0 variable
    %         zp_delta = real(A_c(floor((2*MM+1)/2)+delta,floor((2*NN+1)/2)));
    %         xp_delta = real(A_c(floor((2*MM+1)/2),floor((2*NN+1)/2)+delta));
    %         ZX_0 = real(A_c(floor((2*MM+1)/2),floor((2*NN+1)/2)));
    %         % Estimate shear wave speed for z and x directions and average them
    %         kz(kk,ll) = sqrt(abs(1-zp_delta*T_0/(ZX_0*T_delta_z))*(c_z/delta_z^2));
    %         kx(kk,ll) = sqrt(abs(1-xp_delta*T_0/(ZX_0*T_delta_x))*(c_x/delta_x^2));
    %     end
    % end
pv_field = M_filtered1;  

mirror_pv = padarray(pv_field, (win-1)/2, 'symmetric','both');

step = 2;
og_size = size(M_filtered1);
kz = zeros(length(1:step:og_size(1)),length(1:step:og_size(2)));
kx = zeros(length(1:step:og_size(1)),length(1:step:og_size(2)));

    search_z = 1:step:og_size(1); % To save time, the windowed matrix is evaluated in points separated by a step of 2
    search_x = 1:step:og_size(2);

    for i = 1:length(search_z)
        search_x = 1:step:og_size(2);
        % for j = 1:length(ev_al_pos_x)
        parfor j = 1:length(search_x)
            % Extract relevant submatrix from vz for current position
            % sub_vz = vz(ev_al_pos_z(i) + search_area_z, ev_al_pos_x(j) + search_area_x);
            
            sub_vz = mirror_pv(search_z(i):search_z(i)+win(1)-1,search_x(j):search_x(j)+win(2)-1); % Window kernel EMZ style
            A_c = xcorr2(sub_vz);

             % Define zp and xp vectors and since z(0)= x(0), then they are defined as
            % ZX_0 variable
            zp_delta = real(A_c(floor((2*MM+1)/2)+delta,floor((2*NN+1)/2)));
            xp_delta = real(A_c(floor((2*MM+1)/2),floor((2*NN+1)/2)+delta));
            ZX_0 = real(A_c(floor((2*MM+1)/2),floor((2*NN+1)/2)));
            % Estimate shear wave speed for z and x directions and average them
            kz(i,j) = sqrt(abs(1-zp_delta*T_0/(ZX_0*T_delta_z))*(c_z/delta_z^2));
            kx(i,j) = sqrt(abs(1-xp_delta*T_0/(ZX_0*T_delta_x))*(c_x/delta_x^2));
        end
    end

    k_mean = (kz+kx)/2;
    cs = 2*pi*freq./k_mean;

% Apply a filter to the shear wave speed estimates
mask0 = fspecial('average', [5 5]);
SWS2 = imfilter(cs, mask0, 'replicate');
figure,imagesc(1e2*x_v,1e2*z_v,SWS2,[0 5.2]),colorbar, colormap jet;axis equal tight


sws_big = bigImg(cs,Bmode);
[hF,hB,hColor] = imOverlayInterp(Bmode,sws_big,visDefault.caxis_bmode, visDefault.caxis_img, ...
                0.5,xdim*1e3,ydim*1e3,sws_big,xdim*1e3,ydim*1e3);


%%
%% Estimaging wave velocity

window_size = 51; % must be an odd number
win = [window_size, window_size];
u1=JO_filtered;

mirror_pv = padarray(u1, (win-1)/2, 'symmetric','both');


dz=4; delta_z=dinf.dz*dz;
dx=3; delta_x=dinf.dx*dx;
d135=4; delta_135=dinf.dz*d135;
d45=4; delta_45=dinf.dz*d45;
A_n=xcorr2(ones(window_size,window_size)); % For normalizing autocorrelation

kz = zeros(length(1:step:og_size(1)),length(1:step:og_size(2)));
kx = zeros(length(1:step:og_size(1)),length(1:step:og_size(2)));
sws_fer = kz;

search_z = 1:step:og_size(1); % To save time, the windowed matrix is evaluated in points separated by a step of 2
search_x = 1:step:og_size(2);

    for i = 1:length(search_z)
        search_x = 1:step:og_size(2);
        % for j = 1:length(ev_al_pos_x)
        for j = 1:length(search_x)

        % work_area=z_pad(b:window_size+b-1,a:window_size+a-1);
        work_area = mirror_pv(search_z(i):search_z(i)+win(1)-1,search_x(j):search_x(j)+win(2)-1); % Window kernel EMZ style
       
        auto_corr=xcorr2(work_area);
        auto_corr=auto_corr./A_n;
        
        z_vector=transpose(auto_corr(:,window_size));
        z_vector=z_vector./max(abs(z_vector)); % Normalizing
        k_z2=(5/((delta_z^2)*real(z_vector(window_size))))*(real(z_vector(window_size))-real(z_vector(window_size-dz)));
        
        x_vector=auto_corr(window_size,:);
        x_vector=x_vector./max(abs(x_vector)); % Normalizing
        k_x2=(5/(sqrt(2*pi)*(delta_x^2)*real(x_vector(window_size))))*(real(x_vector(window_size))-real(x_vector(window_size-dx)));
        
        vector_45=transpose(diag(fliplr(auto_corr)));
        vector_45=vector_45./max(abs(vector_45)); % Normalizing
        k_45=(5/(sqrt(2*pi)*(delta_45^2)*real(vector_45(window_size))))*(real(vector_45(window_size))-real(vector_45(window_size-d45)));
        
        vector_135=transpose(diag(auto_corr));
        vector_135=vector_135./max(abs(vector_135)); % Normalizing
        k_135=(5/(sqrt(2*pi)*(delta_135^2)*real(vector_135(window_size))))*(real(vector_135(window_size))-real(vector_135(window_size-d135)));
        
        k_final=sqrt((k_z2+k_x2+k_45+k_135)/4);
        if a==fix(size(u,2)/2) && b==fix(size(u,1)/2)
            figure
            plot(dinf.dz*(10^3)*(-(window_size-1):window_size-1),real(z_vector)); hold on;
            plot(dinf.dx*(10^3)*(-(window_size-1):window_size-1),real(x_vector)); hold off;
            grid on
            xlabel('size [mm]');ylabel('Normalized autocorrelation function');legend('axial','lateral'); title('Autocorrelation profiles');
            figure
            x_scale_au=linspace(-window_size,window_size)*dinf.dx*(10^3); % in mm
            z_scale_au=linspace(-window_size,window_size)*dinf.dz*(10^3); % in mm
            imagesc(x_scale_au,z_scale_au,abs(auto_corr)); colormap jet; colorbar; 
            xlabel('Lateral [mm]');ylabel('Axial [mm]');title('Autocorrelation Map'); 
        end
        sws_fer(i,j)=(2*pi*freq)/k_final;
    end
end


sws_big = bigImg(sws_fer,Bmode);
[hF,hB,hColor] = imOverlayInterp(Bmode,sws_big,visDefault.caxis_bmode, visDefault.caxis_img, ...
                0.5,xdim*1e3,ydim*1e3,sws_big,xdim*1e3,ydim*1e3);