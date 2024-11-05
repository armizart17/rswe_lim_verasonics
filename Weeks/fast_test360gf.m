%%% DATA SIMUS

% C:\Users\emirandaz\OneDrive - pucp.pe\Documentos\MATLAB\TESIS_MPSID\R-SWE\phasegradientRSWE\dataold\Data1000Hz-1
load('C:\Users\emirandaz\OneDrive - pucp.pe\Documentos\MATLAB\TESIS_MPSID\R-SWE\phasegradientRSWE\dataold\Data1000Hz-10000wvs\R-FIELD_inc_1.mat')

% load('C:\Users\emirandaz\OneDrive - pucp.pe\Documentos\MATLAB\TESIS_MPSID\R-SWE\phasegradientRSWE\dataold\Data500Hz-10000wvs\R-FIELD_inc_1.mat')

% load('C:\Users\emirandaz\OneDrive - pucp.pe\Documentos\MATLAB\TESIS_MPSID\R-SWE\phasegradientRSWE\dataold\Data700Hz-10000wvs\R-FIELD_inc_1.mat')

%%

clc;
window = 15; %11 pixels as described in paper
w_kernel = [window window];
stride = 2;

dinf.dx = min(diff(x));
dinf.dz = min(diff(z));
frame = pv_complexZ(:,:,1); % number of frame
        
%frame = (frame'); %transpose for Z (vertical-axial) X(horizontal-lateral)
        
og_size = size(frame);
mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');

% FUNCION SWS MAP (Con linearizacion)
tic
[grad_abs, size_out] = pg_norm(mirror_frame, w_kernel, dinf, og_size, stride);
t = toc;
fprintf("Time PG: %.2f secs \n", t)

tic
num_angles = 3;
[grad_abs_360, size_out_360] = pg_norm360(mirror_frame, w_kernel, dinf, og_size, stride, num_angles);
t = toc;
fprintf("Time PG-360º (%d-angles): %.2f secs \n", num_angles, t)

%%
figure, 
subplot(1,2,1)
imagesc(x*1e3,z*1e3, grad_abs)
title('PG')

subplot(1,2,2)
imagesc(x*1e3,z*1e3, grad_abs_360)
title(['PG-360 ', num2str(num_angles), '-angles'])


% subplot(133)
% imagesc(x*1e3,z*1e3, grad_avg_360)
% title(['PG-360 ', num2str(num_angles), '-angles'])


%
avg_kernel = ones(9, 9) / 81;  % Create a 7x7 averaging filter kernel

grad_abs_avg = filter2(avg_kernel, grad_abs, 'same');
sws_pg = 2*pi*freq ./ grad_abs_avg;

grad_abs_360_avg = filter2(avg_kernel, grad_abs_360, 'same');
sws_pg_360 = 2*pi*freq ./ grad_abs_360_avg;


%
figure, 
subplot(121)
imagesc(x*1e3,z*1e3, sws_pg, [0 6]), colormap("jet")
colorbar
axis("tight")
title('PG')
xlabel('Lateral [1e3]'), ylabel('Depth [1e3]')

subplot(122)
imagesc(x*1e3,z*1e3, sws_pg_360, [0 6]), colormap("jet")
title(['PG-360 ', num2str(num_angles), '-angles'])
colorbar
axis("tight")
xlabel('Lateral [1e3]'), ylabel('Depth [1e3]')

%%
tic
num_angles = 6;
[grad_abs_360_6, size_out_360] = pg_norm360(mirror_frame, w_kernel, dinf, og_size, stride, num_angles);
t = toc;
fprintf("Time PG-360º (%d-angles): %.2f secs \n", num_angles, t)

tic
num_angles = 12;
[grad_abs_360_12, size_out_360] = pg_norm360(mirror_frame, w_kernel, dinf, og_size, stride, num_angles);
t = toc;
fprintf("Time PG-360º (%d-angles): %.2f secs \n", num_angles, t)

tic
num_angles = 24;
[grad_abs_360_24, size_out_360] = pg_norm360(mirror_frame, w_kernel, dinf, og_size, stride, num_angles);
t = toc;
fprintf("Time PG-360º (%d-angles): %.2f secs \n", num_angles, t)

tic
num_angles = 36;
[grad_abs_360_36, size_out_360] = pg_norm360(mirror_frame, w_kernel, dinf, og_size, stride, num_angles);
t = toc;
fprintf("Time PG-360º (%d-angles): %.2f secs \n", num_angles, t)

%%

grad_abs_360_avg = filter2(avg_kernel, grad_abs_360_6, 'same');
sws_pg_360_6 = 2*pi*freq ./ grad_abs_360_avg;

grad_abs_360_avg = filter2(avg_kernel, grad_abs_360_12, 'same');
sws_pg_360_12 = 2*pi*freq ./ grad_abs_360_avg;

grad_abs_360_avg = filter2(avg_kernel, grad_abs_360_24, 'same');
sws_pg_360_24 = 2*pi*freq ./ grad_abs_360_avg;

grad_abs_360_avg = filter2(avg_kernel, grad_abs_360_36, 'same');
sws_pg_360_36 = 2*pi*freq ./ grad_abs_360_avg;

figure, 

subplot(2,2,1),
num_angles = 6;
imagesc(x*1e3,z*1e3, sws_pg_360_6, [0 6]), colormap("jet")
title(['PG-360 ', num2str(num_angles), '-angles'])
colorbar
axis("tight")
xlabel('Lateral [1e3]'), ylabel('Depth [1e3]')

subplot(2,2,2),
num_angles = 12;
imagesc(x*1e3,z*1e3, sws_pg_360_12, [0 6]), colormap("jet")
title(['PG-360 ', num2str(num_angles), '-angles'])
colorbar
axis("tight")
xlabel('Lateral [1e3]'), ylabel('Depth [1e3]')

subplot(2,2,3),
num_angles = 24;
imagesc(x*1e3,z*1e3, sws_pg_360_24, [0 6]), colormap("jet")
title(['PG-360 ', num2str(num_angles), '-angles'])
colorbar
axis("tight")
xlabel('Lateral [1e3]'), ylabel('Depth [1e3]')

subplot(2,2,4),
num_angles = 36;
imagesc(x*1e3,z*1e3, sws_pg_360_36, [0 6]), colormap("jet")
title(['PG-360 ', num2str(num_angles), '-angles'])
colorbar
axis("tight")
xlabel('Lateral [1e3]'), ylabel('Depth [1e3]')


%% 
% PG SIMPLE ROTATION 21/10/2024

clc;
window = 15; %11 pixels as described in paper
w_kernel = [window window];
stride = 2;

dinf.dx = min(diff(x));
dinf.dz = min(diff(z));
pv_field = pv_complexZ(:,:,1); % number of frame
              
og_size = size(pv_field);
pv_field_pad = padarray(pv_field, (w_kernel-1)/2, 'symmetric');

% PG : LINE 171
% PG-ROTATION: LINE 410
%% PG SPLIT
stride = 2;
u = pv_field_pad;

    st = stride;
    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resolution
    
     % Axis of the kernels
    z_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(1))*res_z; 
    x_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(2))*res_x; 
    
    angle_u = angle(u);
    
    size_mirror = size(u); % ogsize + w_kernel - 1; %  = size (u)
    numkernels = floor( (size_mirror - w_kernel)./st + 1 ); % numKernels
    size_out = floor( (og_size - 1)./st + 1 );
   	     
    cont_kernel = 1; 
    angle_z = unwrap(angle_u,[],1);
    angle_x = unwrap(angle_u,[],2);

    angle_2d = unwrap_phase(angle_u); % v2

    grad_l2 = zeros(size_out);
    
    % NORMAL LOOP
    i_out = 1;
    for ii = 1:st:og_size(1)

        j_out = 1;

        for jj = 1:st:og_size(2) %% for faster computing pararell toolbox

            % V1 UNWRAPP XZ
            % area_z = angle_z(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            % area_x = angle_x(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel

            % V2 UNWRAPP 2D
            % area_z = angle_2d(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            % area_x = angle_2d(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel

            % V3 UNWRAPP XZ PER KERNEL
            area_z = angle_u(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            area_x = angle_u(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel

            area_z = unwrap(area_z, [], 1);
            area_x = unwrap(area_x, [], 2);

            % V4 UNWRAPP 2D PER KERNEL
            area_z = angle_u(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            % area_x = angle_u(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel

            area_z = unwrap_phase(area_z);
            % area_x = unwrap_phase(area_x); 
            area_x = area_z; % faster

            diff_area_z = diff(area_z, 1, 1) ./ res_z;
            diff_area_x = diff(area_x, 1, 2) ./ res_x;

            diff_area_z = [diff_area_z; diff_area_z(end, :)];  % replicate the last row
            diff_area_x = [diff_area_x, diff_area_x(:, end)];  % replicate the last column

            val_abs = sqrt(diff_area_z.^2 + diff_area_x.^2);
            grad_l2(i_out,j_out) = mean(val_abs, 'all');

            j_out = j_out + 1;

        end
       i_out = i_out + 1;
    end
      
%% AVG KERNEL
avg_kernel = ones(9, 9) / 81;  % Create a 7x7 averaging filter kernel

caxis_sws = [0 5];

%% V1 UNWRAPP XZ 
grad_l2_avg = filter2(avg_kernel, grad_l2, 'same');

unwrappXZ.grad_l2 = grad_l2; % Unwrapp XZ

unwrappXZ.sws = 2*pi*freq ./ grad_l2_avg;

%% V2 UNWRAPP 2D
grad_l2_avg = filter2(avg_kernel, grad_l2, 'same');

unwrapp2D.grad_l2 = grad_l2; % Unwrapp 2D

unwrapp2D.sws = 2*pi*freq ./ grad_l2_avg;

%% V3 UNWRAPP XZ PER KERNEL

grad_l2_avg = filter2(avg_kernel, grad_l2, 'same');

unwrappXZkernel.grad_l2 = grad_l2; % Unwrapp 2D

unwrappXZkernel.sws = 2*pi*freq ./ grad_l2_avg;

%% V4 UNWRAPP 2D PER KERNEL

grad_l2_avg = filter2(avg_kernel, grad_l2, 'same');

unwrapp2Dkernel.grad_l2 = grad_l2; % Unwrapp 2D

unwrapp2Dkernel.sws = 2*pi*freq ./ grad_l2_avg;
%% FIGURE COMPARISON

figure, 
tiledlayout(1,4)

nexttile
imagesc(unwrappXZ.sws, caxis_sws), colorbar, colormap("turbo")
axis('image');
title('UnwrappFull XZ'), 

% nexttile
% imagesc(unwrapp2D.sws, caxis_sws), colorbar, colormap("turbo")
% axis('image');
% title('UnwrappFull 2D '),

nexttile
imagesc(unwrappXZkernel.sws, caxis_sws), colorbar, colormap("turbo")
axis('image');
title('UnwrappKernel 2D'),

nexttile
imagesc(unwrapp2Dkernel.sws, caxis_sws), colorbar, colormap("turbo")
axis('image');
title('UnwrappKernel 2D'),

%% METRICS COMPARISON

% BIGGER
unwrappXZ.sws = bigImg(unwrappXZ.sws, pv_field);
unwrapp2D.sws = bigImg(unwrapp2D.sws, pv_field);
unwrappXZkernel.sws = bigImg(unwrappXZkernel.sws, pv_field);
unwrapp2Dkernel.sws = bigImg(unwrapp2Dkernel.sws, pv_field);

cx = 0E-3; cz = 0E-3;
radious = 10E-3;

[XX,ZZ] = meshgrid(x,z);

% '1' inside, '0', outside
circ_mask = sqrt( (XX - cx).^2 + (ZZ - cz).^2 ) <= radious; 


unwrappXZ.sws_in_m = mean(unwrappXZ.sws(circ_mask));
unwrappXZ.sws_in_s = std(unwrappXZ.sws(circ_mask));
unwrappXZ.sws_bg_m = mean(unwrappXZ.sws(~circ_mask));
unwrappXZ.sws_bg_s = std(unwrappXZ.sws(~circ_mask));

unwrapp2D.sws_in_m = mean(unwrapp2D.sws(circ_mask));
unwrapp2D.sws_in_s = std(unwrapp2D.sws(circ_mask));
unwrapp2D.sws_bg_m = mean(unwrapp2D.sws(~circ_mask));
unwrapp2D.sws_bg_s = std(unwrapp2D.sws(~circ_mask));

unwrappXZkernel.sws_in_m = mean(unwrappXZkernel.sws(circ_mask));
unwrappXZkernel.sws_in_s = std(unwrappXZkernel.sws(circ_mask));
unwrappXZkernel.sws_bg_m = mean(unwrappXZkernel.sws(~circ_mask));
unwrappXZkernel.sws_bg_s = std(unwrappXZkernel.sws(~circ_mask));

unwrapp2Dkernel.sws_in_m = mean(unwrapp2Dkernel.sws(circ_mask));
unwrapp2Dkernel.sws_in_s = std(unwrapp2Dkernel.sws(circ_mask));
unwrapp2Dkernel.sws_bg_m = mean(unwrapp2Dkernel.sws(~circ_mask));
unwrapp2Dkernel.sws_bg_s = std(unwrapp2Dkernel.sws(~circ_mask));


unwrappXZ.method = 'UnwrappFull XZ';
unwrapp2D.method = 'UnwrappFull 2D';
unwrappXZkernel.method = 'UnwrappKernel XZ';
unwrapp2Dkernel.method = 'UnwrappKernel 2D';

pmSymbol = char(177); %+/-

unwrappXZ.title = string( sprintf('%s, in: %.2f %c %.2f, bg: %.2f %c %.2f' ...
                        , unwrappXZ.method ...
                        , unwrappXZ.sws_in_m, pmSymbol,unwrappXZ.sws_in_s ...
                        , unwrappXZ.sws_bg_m, pmSymbol,unwrappXZ.sws_bg_s ...
                        ) );

unwrapp2D.title = string( sprintf('%s, in: %.2f %c %.2f, bg: %.2f %c %.2f' ...
                        , unwrapp2D.method ...
                        , unwrapp2D.sws_in_m, pmSymbol,unwrapp2D.sws_in_s ...
                        , unwrapp2D.sws_bg_m, pmSymbol,unwrapp2D.sws_bg_s ...
                        ) );

unwrappXZkernel.title = string (sprintf('%s, in: %.2f %c %.2f, bg: %.2f %c %.2f' ...
                        , unwrappXZkernel.method ...
                        , unwrappXZkernel.sws_in_m, pmSymbol,unwrappXZkernel.sws_in_s ...
                        , unwrappXZkernel.sws_bg_m, pmSymbol,unwrappXZkernel.sws_bg_s ...
                        ) );

unwrapp2Dkernel.title = sprintf('%s, in: %.2f %c %.2f, bg: %.2f %c %.2f' ...
                        , unwrapp2Dkernel.method ...
                        , unwrapp2Dkernel.sws_in_m, pmSymbol,unwrapp2Dkernel.sws_in_s ...
                        , unwrapp2Dkernel.sws_bg_m, pmSymbol,unwrapp2Dkernel.sws_bg_s ...
                        ) ;

%%
figure, 
tiledlayout(1,4)
sgtitle(['NumAngles=', num2str(num_angles)]); % FOR ANGLES

nexttile
imagesc(x*1e3, z*1e3, unwrappXZ.sws, caxis_sws), colorbar, colormap("turbo")
axis('image');
title(unwrappXZ.title)

nexttile
imagesc(x*1e3, z*1e3, unwrapp2D.sws, caxis_sws), colorbar, colormap("turbo")
axis('image');
title(unwrapp2D.title)

nexttile
imagesc(x*1e3, z*1e3, unwrappXZkernel.sws, caxis_sws), colorbar, colormap("turbo")
axis('image');
title(unwrappXZkernel.title)

nexttile
imagesc(x*1e3, z*1e3,unwrapp2Dkernel.sws, caxis_sws), colorbar, colormap("turbo")
axis('image');
title(unwrapp2Dkernel.title)

%%

T = [
    struct2table(unwrappXZ);
    struct2table(unwrapp2D);
    struct2table(unwrappXZkernel);
    struct2table(unwrapp2Dkernel)
    ];

%%
writetable(T,fullfile(pathout,'results.xlsx'),'WriteRowNames',true);
close all


%% NOW PG-ROTATION ATTEMPT
stride = 2;
u = pv_field_pad;
num_angles = 12;

    tic;

    % Inputs parsing and initial setup
    st = stride;
    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resolution

    % Axis of the kernels
    z_axis = linspace(-(w_kernel(1)-1)/2, (w_kernel(1)-1)/2, w_kernel(1)) * res_z;
    x_axis = linspace(-(w_kernel(2)-1)/2, (w_kernel(2)-1)/2, w_kernel(2)) * res_x;
    
    % Initializing phase angles and output size
    angle_u = angle(u);
    size_mirror = size(u);
    numkernels = floor((size_mirror - w_kernel) ./ st + 1);
    size_out = floor((og_size - 1) ./ st + 1);

    angle_z = unwrap(angle_u, [], 1);
    angle_x = unwrap(angle_u, [], 2);
    
    angle_2d = unwrap_phase(angle_u); % v2

    % figure, 
    % subplot(2,2,1), imagesc(angle_u); title('angle u')
    % subplot(2,2,3), imagesc(angle_z); title('angle z')
    % subplot(2,2,4), imagesc(angle_x); title('angle x')
    % subplot(2,2,2), imagesc(angle_2d); title('angle 2d')

    grad_l2 = zeros(size_out);
    grad_l2_circ = zeros(size_out);

    % Compute the search space once
    search_z = 1:st:og_size(1);
    search_x = 1:st:og_size(2);

    % Calculate the rotation angle step based on the number of angles
    delta_theta = 360 / num_angles;  % This might be a decimal value

    % Pre-allocate 3D array for storing gradients across rotations
    val_abs_3D = zeros(w_kernel(1), w_kernel(2), num_angles);

    % Circle mask for rotation 
    M = w_kernel(1); N = w_kernel(2);
    [xGrid, yGrid] = meshgrid(1:N, 1:M);

    R = floor (min(M, N) / 2) ;
    centerX = ceil(N / 2);
    centerY = ceil(M / 2);

    circ_mask  = sqrt((xGrid - centerX).^2 + (yGrid - centerY).^2) < R;


    % Main loop for gradient calculation
    for ii = 1:length(search_z)
        for jj = 1:length(search_x)
            
            % V1 UNWRAPP XZ
            % area_z = angle_z(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1); 
            % area_x = angle_x(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1);

            % V2 UNWRAPP 2D
            % area_z = angle_2d(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1); 
            % area_x = area_z;

            % V3 UNWRAPP XZ PER KERNEL      
            % area_z = angle_u(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1); 
            % area_x = angle_u(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1);
            % 
            % area_z = unwrap(area_z, [], 1);
            % area_x = unwrap(area_x, [], 2);

            % V4 UNWRAPP 2D PER KERNEL
            area_2d = angle_u(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1); 

            area_z = unwrap_phase(area_2d);
            area_x = area_z;

            % if(ii == fix( length(search_z) / 2) && ...
            % disp(ii)
            % if(search_z(ii) == 271 && ... % equivale a 4cm depth             
            % search_x(jj) == ceil( length(search_x) / 2) )
            %     keyboard;
            % end
            % 
            % if(search_z(ii) == 271 && ... % equivale a 4cm depth             
            % search_x(jj) == ceil( length(search_x) / 3) )
            %     keyboard;
            % end

            if(search_z(ii) == 271 && ... % equivale a 4cm depth             
            search_x(jj) == 21 ) %-1.28cm
                keyboard;
            end

            if(search_z(ii) == 271 && ... % equivale a 4cm depth             
            search_x(jj) == 101 ) % 1.08 cm 
                keyboard;
            end

            % Parallel loop over rotations using the given number of angles
            for idx = 1:num_angles

                clear diff_area_z diff_area_x
                % Calculate the actual rotation angle
                theta = mod((idx - 1) * delta_theta, 360); % Ensuring the angle wraps around at 360

                % Rotate areas and calculate gradients
                rotated_area_z = imrotate(area_z, theta, 'nearest', 'crop');
                rotated_area_x = imrotate(area_x, theta, 'nearest', 'crop');
                % 
                % figure, 
                % subplot(221), imagesc(area_x), title('x'), colorbar
                % subplot(222), imagesc(area_z), title('z'), colorbar
                % subplot(223), imagesc(rotated_area_x), title('Rot x'), colorbar
                % subplot(224), imagesc(rotated_area_z), title('Rot z'), colorbar

                diff_area_z = diff(rotated_area_z, 1, 1) ./ res_z;
                diff_area_x = diff(rotated_area_x, 1, 2) ./ res_x;

                % Ensure consistent dimensions by padding the last row/column
                diff_area_z = [diff_area_z; diff_area_z(end, :)]; % replicate last row
                diff_area_x = [diff_area_x, diff_area_x(:, end)];  % replicate last column

                % Compute the L2-norm of the gradients
                val_abs_3D(:, :, idx) = sqrt(diff_area_z.^2 + diff_area_x.^2);

            end

            % Average L2-norm across all rotations through depth (angle dimension)
            val = mean(val_abs_3D, 3); % 2D-array

            val_circ_mask = val(circ_mask); % select circle mask returns a 1D-array


            % grad_l2(ii, jj) = mean(val, 'all');
            grad_l2(ii,jj) = mean(val_circ_mask,'all');

        end
    end

    tt = toc;
    fprintf('Execution time PG-Rot %.4f \n', tt);

%% AVG KERNEL
avg_kernel = ones(9, 9) / 81;  % Create a 7x7 averaging filter kernel

caxis_sws = [0 5];

%% V1 UNWRAPP XZ 
grad_l2_avg = filter2(avg_kernel, grad_l2, 'same');

unwrappXZ.grad_l2 = grad_l2; % Unwrapp XZ

unwrappXZ.sws = 2*pi*freq ./ grad_l2_avg;

%% V2 UNWRAPP 2D
grad_l2_avg = filter2(avg_kernel, grad_l2, 'same');

unwrapp2D.grad_l2 = grad_l2; % Unwrapp 2D

unwrapp2D.sws = 2*pi*freq ./ grad_l2_avg;

%% V3 UNWRAPP XZ PER KERNEL

grad_l2_avg = filter2(avg_kernel, grad_l2, 'same');

unwrappXZkernel.grad_l2 = grad_l2; % Unwrapp 2D

unwrappXZkernel.sws = 2*pi*freq ./ grad_l2_avg;

%% V4 UNWRAPP 2D PER KERNEL

grad_l2_avg = filter2(avg_kernel, grad_l2, 'same');

unwrapp2Dkernel.grad_l2 = grad_l2; % Unwrapp 2D

unwrapp2Dkernel.sws = 2*pi*freq ./ grad_l2_avg;