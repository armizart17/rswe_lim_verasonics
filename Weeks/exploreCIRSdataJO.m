% simple code data JO data2_PUCP 
addpath(genpath(pwd));
fprintf('------PG l2 norm -------\n')

list_data = [10 11 13 14];
v_freq = [400 600 900 900];

stride_pg = 2;
w_kernel_pg = [15 15];
tic;

pathdata = 'P:/rswe/Data_3_PUCP/CIRS_phantom/L7-4/';

for ii = 1: length(v_freq)
   
    freq = v_freq(ii);
    sampleCase = list_data(ii);
    pathfreq_in = [pathdata,'data_', num2str(sampleCase),'/'];

    name = ['data_',num2str(sampleCase),'.mat'];
    structdata = load([pathfreq_in, name]);

    dinf.dx = structdata.dinf.dx;
    dinf.dz = structdata.dinf.dz;
    dinf.PRFe = structdata.dinf.PRFe;

    % NORMALY BUT WE ALREADY HAVE IT
    % [u,dinf] = pv_cal(IQData1,dinf, dinf.num_angles); 
    % u = signal_period(freq, dinf.PRFe, u);  

    % 
    [JO_filtered] = fun_JO_v1(structdata.u, freq, dinf);

    pv = JO_filtered;
    og_size = size(pv);
    mirror_pv = padarray(pv,(w_kernel_pg-1)/2,'symmetric');

    % FUNCION PG-l2 norm
    [grad_l2, size_out] = pg_norm(mirror_pv, w_kernel_pg, dinf, og_size, stride_pg);
    
    pg.grad_l2_3D (:, :, ii) = grad_l2;
    pg.sws_l2_3D(:,:, ii) = 2*pi*freq ./ grad_l2;
    pg.list_freq(ii) = freq;

end

%%
%% PG DENOISING WITH AVERAGE FILTER
sws_max = 5; % [m/s]
caxis_sws = [0 sws_max];
numChannels = length(pg.list_freq);
nCols = ceil (numChannels/2);

size_avg = 7;
Bmode = 20*log10(abs(structdata.IQBmodeData(13:end,:,1)));
Bmode = Bmode - max(Bmode(:));

xdim = linspace(-dinf.dx*size(Bmode,2)/2,dinf.dx*size(Bmode,2)/2,size(Bmode,2)); 
ydim = linspace(0,dinf.dz*size(Bmode,1),size(Bmode,1));

x = xdim; z = ydim;

for ii = 1 : numChannels
    freq = pg.list_freq(ii);
    avg_kernel = ones(size_avg, size_avg) / size_avg^2;  % Create a 7x7 averaging filter kernel
    avef.grad_l2_3D(:,:,ii) = filter2(avg_kernel, pg.grad_l2_3D(:,:,ii), 'same');
    avef.sws_l2_3D(:,:,ii) = 2*pi*freq ./ avef.grad_l2_3D(:, :, ii);
end

% PLOT AVERAGE FILT 

figure, % SWS
sgtitle('\bfSWS Average filter')
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.75])
for ii = 1 : numChannels
    freq = pg.list_freq(ii);
    subplot (2, nCols, ii)
    imagesc(x*1e2, z*1e2, avef.sws_l2_3D(:,:,ii), caxis_sws)
%         colormap ("jet");
    colormap ("turbo");
    % axis("equal")
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS @ ', num2str(freq), 'Hz' ])

end
%%
% simple code data JO data2_PUCP 
addpath(genpath(pwd));
fprintf('------CF -------\n')

list_data = [10 11 13 14];
v_freq = [400 600 900 900];

stride_cf = 2;
w_kernel_cf = [101 101];
tic;

pathdata = 'P:/rswe/Data_3_PUCP/CIRS_phantom/L7-4/';

for ii = 1: length(v_freq)
   
    freq = v_freq(ii);
    sampleCase = list_data(ii);
    pathfreq_in = [pathdata,'data_', num2str(sampleCase),'/'];


    name = ['data_',num2str(sampleCase),'.mat'];
    structdata = load([pathfreq_in, name]);

    dinf.dx = structdata.dinf.dx;
    dinf.dz = structdata.dinf.dz;

    % NORMALY BUT WE ALREADY HAVE IT
    % [u,dinf] = pv_cal(IQData1,dinf, dinf.num_angles); 
    % u = signal_period(freq, dinf.PRFe, u);  

    % 
    [JO_filtered] = fun_JO_v1(structdata.u, freq, dinf);

    pv = JO_filtered;
    og_size = size(pv);
    mirror_pv = padarray(pv,(w_kernel_cf-1)/2,'symmetric');

    % FUNCION CF

    correc = xcorr2(ones(w_kernel_cf(1), w_kernel_cf(2))); % it is default 
    tic 
    [Kx,Kz,Rx,Rz,K1d,R1d] = sws_estimation_cf_fast(mirror_pv, w_kernel_cf, dinf.dx, dinf.dz, correc, og_size); % EMZ
    tt = toc;
    fprintf('Time passed CF %.4f\n', tt)
    
    cf.Kx (:, :, ii) = Kx;
    cf.Kz (:, :, ii) = Kz;
    cf.Rx (:, :, ii) = Rx;
    cf.Rz (:, :, ii) = Rz;
    cf.K1d (:, :, ii) = K1d;
    cf.R1d (:, :, ii) = R1d;
    cf.list_freq(ii) = freq;

end

%%
%%
% simple code data JO data2_PUCP 
addpath(genpath(pwd));
fprintf('------MAOW -------\n')

list_data = [10 11 13 14];
v_freq = [400 600 900 900];

stride_cf = 2;
w_kernel_maow = [101 101];
tic;

pathdata = 'P:/rswe/Data_3_PUCP/CIRS_phantom/L7-4/';

for ii = 1: length(v_freq)
   
    freq = v_freq(ii);
    sampleCase = list_data(ii);
    pathfreq_in = [pathdata,'data_', num2str(sampleCase),'/'];


    name = ['data_',num2str(sampleCase),'.mat'];
    structdata = load([pathfreq_in, name]);

    dinf.dx = structdata.dinf.dx;
    dinf.dz = structdata.dinf.dz;

    % NORMALY BUT WE ALREADY HAVE IT
    % [u,dinf] = pv_cal(IQData1,dinf, dinf.num_angles); 
    % u = signal_period(freq, dinf.PRFe, u);  

    % 
    [JO_filtered] = fun_JO_v1(structdata.u, freq, dinf);

    pv = JO_filtered;
    og_size = size(pv);
    mirror_pv = padarray(pv,(w_kernel_maow-1)/2,'symmetric');

    % FUNCION MAOW
    tic 
    [sws_matrix,k_z,k_x] = sws_generator(u,w_kernel_maow,freq,d_n,dinf,og_size,10,5);
    tt = toc;
    fprintf('Time passed MAOW %.4f\n', tt)

    maow.sws (:, :, ii) = sws_matrix;
    maow.Kz (:, :, ii) = k_z;
    maow.Rx (:, :, ii) = k_x;


end


%%  PG DENOISING WITH TOTAL VARIATION

M = size_out(1);
N = size_out(2);
x = xdim;
z = ydim;

mu = 10^4.;
for ii = 1 : numChannels
    freq = pg.list_freq(ii);

    my_grad = pg.grad_abs_3D (:, :, ii);

    [grad_tv] = IRLS_TV_simple(my_grad(:),speye(M*N),mu,M,N,1e-4,ones(size(M*N)),ones(M*N,1));

    tv.grad_abs_3D(:,:,ii) = reshape(grad_tv, [ M N ] );

    tv.sws_abs_3D(:,:,ii) = (2*pi*freq)./tv.grad_abs_3D(:,:,ii);

end

% PLOT  TV
figure, % SWS
% caxis_sws = [2 5];

% sgtitle(['SWS TV, \mu=', num2str(mu)])
sgtitle('\bfSWS TV')
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
for ii = 1 : numChannels
    freq = pg.list_freq(ii);
    subplot (2, nCols, ii)
    imagesc(x*1e2, z*1e2, tv.sws_abs_3D(:,:,ii), caxis_sws)
%     colormap ("jet");
    colormap("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS @ ', num2str(freq),'Hz' ])

end