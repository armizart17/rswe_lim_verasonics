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
%% PHASE GRADIENT L2-norm
methodName = 'PG-L2';
stride = 2; % Stride for window-based method
% window = 31;

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

figure, 
imagesc(sws_pg)
title(methodName)

%% PHASE GRADIENT WITH TOTAL VARIATION (PG-TV)
methodName = 'PG-TV';

% Fix to 0.5 wv 

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

figure, 
imagesc(sws_pg_tv)
title(methodName)
