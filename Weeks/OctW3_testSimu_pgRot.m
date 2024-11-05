% Script to test new pg_norm_rot estimator
% Apply a rotation to each kernel
% Simulation case circular inclusion

% (i) Determine effect of angles 
% (ii) Best angle at different frequencies

%% SIMULATION 
pathData = 'C:\Users\emirandaz\OneDrive - pucp.pe\Documentos\MATLAB\TESIS_MPSID\R-SWE\phasegradientRSWE\dataold\';


freq = 500;

load(fullfile(pathData,"Data"+ freq +"Hz-10000wvs\R-FIELD_inc_1.mat"));

window = 15; %11 pixels as described in paper
w_kernel = [window window];
stride = 2;

dinf.dx = min(diff(x));
dinf.dz = min(diff(z));
pv_field = pv_complexZ(:,:,1); % number of frame
              
og_size = size(pv_field);
pv_field_pad = padarray(pv_field, (w_kernel-1)/2, 'symmetric');


%%
clear pg
l_numAngles = [1 3 4 6 9 ];
version = 4;

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
sws_max = 5; % [m/s]
caxis_sws = [0 sws_max];
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

    [grad_tv] = IRLS_TV_simple(my_grad(:),speye(M*N),mu,M,N,1e-4,ones(size(M*N)),ones(M*N,1));

    tv.grad_l2_angle(:,:,ii) = reshape(grad_tv, [ M N ] );

    tv.sws_angle(:,:,ii) = (2*pi*freq)./tv.grad_l2_angle(:,:,ii);

end

% PLOT  TV
figure, % SWS
caxis_sws = [0 5];

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

%% METRICS

%% REGION METRICS CIRCULAR FULL
% 
% IDEAL CONTOURS
radius = 10e-3; centerXZ = [0 0]*1e-2; 
x_circ = centerXZ(1) - radius; w_circ= centerXZ(1) * 2;
y_circ = centerXZ(2) - radius; h_circ = centerXZ(2) * 2;


% EXTRACT MASKS
refBig_sws = pv_complexZ(:,:,1);

% INCLUSION MASK
r_inc = 10E-3; % slightly smaller
[~, inc] = mask_circ(x, z, centerXZ(1), centerXZ(2), r_inc, refBig_sws);

% BACKGROUND MASK
r_bac = 10E-3; % slightly bigger
[~, back] = mask_circ(x, z, centerXZ(1), centerXZ(2), r_bac, refBig_sws);
back = ~back;

% figure, 
% subplot(1,2,1), imagesc(x*1e3, z*1e3,inc), title('Mask Inc'), colorbar, axis('tight');
% subplot(1,2,2), imagesc(x*1e3, z*1e3,back), title('Mask Back'), colorbar,  axis('tight');


gt_inc = 4.5; % [m/s]
gt_back = 2.5; % [m/s]

clear MetricsPG MetricsTV 

% MetricsTV(numChannels) = struct('num_angles', []);
for ii = 1 : numChannels

    sws =  bigImg( tv.sws_angle(:,:,ii) , pv_complexZ(:,:,1) );
    MetricsTV(ii) = get_metrics(sws,inc,back,'tv',freq, gt_inc, gt_back);
   
    % num_angles = l_numAngles(ii);
    % MetricsTV(ii).num_angles = num_angles;
end

% [MetricsTV.num_angles] = deal(l_numAngles);
T = [struct2table(MetricsTV)];

T.num_angles = l_numAngles';

%% Plots DISPERSION INC
% cnr_tv = T.cnr(strcmp(T.method, 'tv'));
lw = 2;
fontSize = 22;
figure,
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.55 0.55])
title('\bfSWS dispersion inclusion')
hold on
errorbar(l_numAngles,T.mean_inc(1:6),T.std_inc(1:6), 'r-o', 'MarkerFaceColor','r','MarkerSize',10, 'LineWidth',lw, 'DisplayName', 'PG-TV')
% errorbar(v_freq,T.mean_inc(7:12),T.std_inc(7:12), 'b-s', 'MarkerFaceColor','b','MarkerSize',10, 'LineWidth',lw)
% errorbar(v_freq,T.mean_inc(13:18),T.std_inc(13:18), 'k-', 'MarkerFaceColor','k','MarkerSize',10, 'Marker','hexagram', 'LineWidth',lw)
yline(4.5, 'k--', 'Ground Truth', 'FontSize',fontSize, 'LineWidth',lw); % 'k--' denotes a black dashed line

hold off
% legend('PG','PG-TV','PG-TNV', 'Location','northeast')
legend('Location', 'northeast')


grid on
% xlim([450 1050])
ylim([0 14])
ylabel('SWS [m/s]'), xlabel('Num Angles')
set(gca, 'FontSize',fontSize)
%% DISPERSION BACK
lw = 1.5;
figure,
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.55 0.55])
title('\bfSWS dispersion background')
hold on

errorbar(v_freq,T.mean_back(1:6),T.std_back(1:6), 'r-o', 'MarkerFaceColor','r','MarkerSize',10, 'LineWidth',lw)
errorbar(v_freq,T.mean_back(7:12),T.std_back(7:12), 'b-s', 'MarkerFaceColor','b','MarkerSize',10, 'LineWidth',lw)
errorbar(v_freq,T.mean_back(13:18),T.std_back(13:18), 'k-', 'MarkerFaceColor','k','MarkerSize',10, 'Marker','hexagram', 'LineWidth',lw)
yline(2.5, 'k--', 'Ground Truth', 'FontSize',fontSize, 'LineWidth',lw); % 'k--' denotes a black dashed line
hold off
legend('PG','PG-TV','PG-TNV', 'Location','northeast')
grid on
xlim([450 1050])

ylabel('SWS [m/s]'), xlabel('Frequency [Hz]')
set(gca, 'FontSize',fontSize)

%%
%% DISPERSION BOTH
colors = [0.8,0.2,0.2; 0.2,0.2,0.8; 0.1 0.7 0.1];

lw = 2;
fontSize = 22;
figure,
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.55 0.55])
title('\bfSWS dispersion')
hold on
errorbar(l_numAngles,T.mean_inc(1:numChannels),T.std_inc(1:numChannels), 'r-o', 'MarkerFaceColor','r', ...
    'Color','r', 'MarkerSize',10, 'LineWidth',lw, 'DisplayName', 'PG-TV')
% errorbar(v_freq,T.mean_inc(7:12),T.std_inc(7:12), 'b-s', 'MarkerFaceColor',colors(2,:), ...
%     'Color',colors(2,:), 'MarkerSize',10, 'LineWidth',lw)
% errorbar(v_freq,T.mean_inc(13:18),T.std_inc(13:18), 'x-', 'MarkerFaceColor',colors(3,:), ...
%     'Color',colors(3,:), 'MarkerSize',10, 'LineWidth',lw)
yline(4.5, 'k--', 'GT_{in}', 'FontSize',fontSize, 'LineWidth',lw); % 'k--' denotes a black dashed line

errorbar(l_numAngles,T.mean_back(1:numChannels),T.std_back(1:numChannels), 'b-s', 'MarkerFaceColor','b', ...
    'Color','b','MarkerSize',10, 'LineWidth',lw, 'DisplayName', 'PG-TV')
% errorbar(v_freq,T.mean_back(7:12),T.std_back(7:12), 'b-s', 'MarkerFaceColor',colors(2,:), ...
%     'Color',colors(2,:),'MarkerSize',10, 'LineWidth',lw)
% errorbar(v_freq,T.mean_back(13:18),T.std_back(13:18), 'x-', 'MarkerFaceColor',colors(3,:), ...
%     'Color',colors(3,:),'MarkerSize',10, 'LineWidth',lw)
yline(2.5, 'k--', 'GT_{bg}', 'FontSize',fontSize, 'LineWidth',lw); % 'k--' denotes a black dashed line

hold off
legend('Location', 'northeast')
grid on
% xlim([450 1050])
ylim([0 8])
ylabel('SWS [m/s]'), xlabel('NumAngles')
set(gca, 'FontSize',fontSize)

%% CNR
figure,
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.55 0.55])
title('\bfCNR')
hold on
plot(l_numAngles,T.cnr(1:numChannels), 'r-o', 'MarkerFaceColor','r', ...
    'Color','r','MarkerSize',10, 'LineWidth',lw, 'DisplayName', 'PG-TV')
% plot(v_freq,T.cnr(7:12), 'b-s', 'MarkerFaceColor','b', ...
%     'Color','b', 'MarkerSize',10, 'LineWidth',lw)
% plot(v_freq,T.cnr(13:18), 'k-', 'MarkerFaceColor','k', ...
%     'Color','k', 'MarkerSize',10, 'LineWidth',lw, 'Marker','hexagram')

hold off
legend('Location','northeast')
grid on
% xlim([450 1050])
ylim([0.5 2])
set(gca, 'YScale', 'log')
set(gca, 'FontSize',fontSize)
xlabel('NumAngles')
% saveas(gcf,fullfile(pathout,'cnrElastic.png'))