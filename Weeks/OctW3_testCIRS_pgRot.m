% Script to test new pg_norm_rot estimator
% Apply a rotation to each kernel
% Simulation case circular inclusion

% (i) Determine effect of angles 
% (ii) Best angle at different frequencies
% Oct 2024
clear all, close all, clc;
%%

pathData = 'C:\Users\emirandaz\OneDrive - pucp.pe\Documentos\MATLAB\TESIS_MPSID\R-SWE\phasegradientRSWE\Data_3_PUCP\CIRS_phantom\L7-4\';

list_data = [10 11 13];
v_freq = [400 600 900 ];

idxCase = 3;

sampleCase = list_data(idxCase);


name = ['data_',num2str(sampleCase),'.mat'];

load(fullfile(pathData, "data_"+ sampleCase,"data_"+ sampleCase +".mat"));
% x = linspace(-dinf.dx*size(Bmode,2)/2,dinf.dx*size(Bmode,2)/2,size(Bmode,2)); 
% z = linspace(0,dinf.dz*size(Bmode,1),size(Bmode,1));
freq = v_freq(idxCase);

window = 15; %11 pixels as described in paper
w_kernel = [window window];
stride = 1;

pv_field = fun_JO_v1(u, freq, dinf);

x = (0:size(pv_field,2)-1)*dinf.dx;
z = (0:size(pv_field,1)-1)*dinf.dz;
x = x - mean(x);
              
og_size = size(pv_field);
pv_field_pad = padarray(pv_field, (w_kernel-1)/2, 'symmetric');



clear pg
l_numAngles = [1 3];
version = 4;

% if iteration angles
for ii = 1:length(l_numAngles)
    num_angles = l_numAngles(ii);
    tic
    [grad_l2, size_out] = pg_norm_rot(pv_field_pad, w_kernel, dinf, og_size, stride, num_angles,version);
    tt = toc;
    fprintf('Exec Time %d-Angles: %.2f secs \n', num_angles, tt)
    pg.grad_l2_angle(:,:,ii) = grad_l2;
    pg.l_numAngles(ii) = num_angles;
end


%% 1 angle (just for check ok)
% tic
% pg_norm(pv_field_pad, w_kernel, dinf, og_size, stride);
pg.grad_l2_angle(:,:,1) = pg_norm(pv_field_pad, w_kernel, dinf, og_size, stride);
% toc

%% PLOT ORIGINAL PG gradient
% sws_max = input('SWS max: ');
sws_max = 4; % [m/s]
caxis_sws = [1 sws_max];
numChannels = length(l_numAngles);

nRows = 2;
nCols = ceil(numChannels/nRows);

% L2-PHASE GRADIENT
figure,
% set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
tiledlayout(nRows, nCols)
sgtitle('\bf L2-Phase Gradient ||\nabla\phi||_2')

for ii = 1 : numChannels
    num_angles = l_numAngles(ii);
    nexttile
    imagesc(x*1e2,z*1e2,pg.grad_l2_angle(:,:,ii))
    axis("image");
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['||\nabla\phi||_2 @ ', num2str(num_angles), '-Angles' ])
    % title (['||\nabla\phi||_2 f_v = ', num2str(freq) ])
end

% SWS
figure,
tiledlayout(nRows, nCols)
sgtitle('\bf PG-og')
% set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])


for ii = 1 : numChannels
    num_angles = l_numAngles(ii);
    sws_aux = (2*pi*freq)./pg.grad_l2_angle(:,:,ii);
    nexttile
    imagesc(x*1e2,z*1e2,sws_aux, caxis_sws)
    % colormap ("jet");
    axis("image");
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS with ', num2str(num_angles),'-Angles' ])
    % title (['SWS f_v = ', num2str(freq) ])

end


%%  PG DENOISING WITH TOTAL VARIATION

M = size_out(1);
N = size_out(2);
numChannels = length(l_numAngles);

mu = 10^4.;
for ii = 1 : numChannels
    
    my_grad = pg.grad_l2_angle(:,:,ii);

    [grad_tv] = IRLS_TV(my_grad(:),speye(M*N),mu,M,N,1e-4,ones(size(M*N)),ones(M*N,1));

    tv.grad_l2_angle(:,:,ii) = reshape(grad_tv, [ M N ] );

    tv.sws_angle(:,:,ii) = (2*pi*freq)./tv.grad_l2_angle(:,:,ii);

end

% PLOT  TV
figure, % SWS
caxis_sws = [1.5 4.5];

nRows = 2;
nCols = ceil(numChannels/nRows);
tiledlayout(nRows, nCols)
% sgtitle(['SWS TV, \mu=', num2str(mu)])
sgtitle('\bfSWS TV')
% set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])


for ii = 1 : numChannels
    num_angles = l_numAngles(ii);
    nexttile
    imagesc(x*1e2, z*1e2, tv.sws_angle(:,:,ii), caxis_sws)
    axis("image");
%     colormap ("jet");
    colormap("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS with ', num2str(num_angles),'-Angles' ])

end



%% ITERATION FREQUENCIES

pathData = 'C:\Users\emirandaz\OneDrive - pucp.pe\Documentos\MATLAB\TESIS_MPSID\R-SWE\phasegradientRSWE\Data_3_PUCP\CIRS_phantom\L7-4\';

list_data = [10 11 13];
v_freq = [400 600 900 ];

% idxCase = 3;
figure

for idxCase = 1:length(v_freq)

sampleCase = list_data(idxCase);


name = ['data_',num2str(sampleCase),'.mat'];

load(fullfile(pathData, "data_"+ sampleCase,"data_"+ sampleCase +".mat"));
% x = linspace(-dinf.dx*size(Bmode,2)/2,dinf.dx*size(Bmode,2)/2,size(Bmode,2)); 
% z = linspace(0,dinf.dz*size(Bmode,1),size(Bmode,1));
freq = v_freq(idxCase);

window = 15; %11 pixels as described in paper
window = 5; % test one time GF
w_kernel = [window window];
stride = 1;

pv_field = fun_JO_v1(u, freq, dinf);

x = (0:size(pv_field,2)-1)*dinf.dx;
z = (0:size(pv_field,1)-1)*dinf.dz;
x = x - mean(x);
              

subplot(1,3,idxCase)
% imagesc(x*1e2, z*1e2, angle(pv_field))
imagesc(angle(pv_field))
% axis ("image")
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
title(['Angle f_v=', num2str(freq)])

og_size = size(pv_field);
pv_field_pad = padarray(pv_field, (w_kernel-1)/2, 'symmetric');


% if iteration angles
% 
% num_angles = 4;
% version = 4;

num_angles = 1; version =1; % ORIGINAL IUS

[grad_l2, size_out] = pg_norm_rot(pv_field_pad, w_kernel, dinf, og_size, stride, num_angles,version);

% fprintf('Exec Time %d-Angles: %.2f secs \n', num_angles, tt)
pg.grad_l2_freq(:,:,idxCase) = grad_l2;
pg.l_freq(idxCase) = freq;

end
%%
M = size_out(1);
N = size_out(2);
numChannels = length(pg.l_freq);

caxis_sws = [1 4];

mu = 10^4.;
for ii = 1 : numChannels
    
    freq = v_freq(ii);

    my_grad = pg.grad_l2_freq(:,:,ii);

    [grad_tv] = IRLS_TV(my_grad(:),speye(M*N),mu,M,N,1e-4,ones(size(M*N)),ones(M*N,1));

    tv.grad_l2_freq(:,:,ii) = reshape(grad_tv, [ M N ] );

    tv.sws_freq(:,:,ii) = (2*pi*freq)./tv.grad_l2_freq(:,:,ii);

end

% PLOT  TV
figure, % SWS
caxis_sws = [1 4];

nRows = 1;
nCols = ceil(numChannels/nRows);
% tiledlayout(nRows, nCols)
% sgtitle(['SWS TV, \mu=', num2str(mu)])
sgtitle('SWS TV')
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])


for ii = 1 : numChannels
    freq = v_freq(ii);
    % nexttile
    subplot (1, numChannels, ii)
    imagesc(x*1e2, z*1e2, tv.sws_freq(:,:,ii), caxis_sws)
    axis("image");
%     colormap ("jet");
    colormap("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ]);

end