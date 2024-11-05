% SIMPLE SCRIPT FOR CLUSTER TEST
% REQUIRES AT LEAST 4GB OF VRAM for NX size 50

clc; clear; close all;

addpath(genpath('/opt/MATLAB Add-Ons/')); 
addpath(genpath(pwd));

disp(pwd);

DATA_CAST = 'gpuArray-single';  % GPU
% DATA_CAST = 'single';           % CPU

iniDir = '/mnt/nfs/emiranda/proj/rswe/'; % for cluster

addpath(genpath(iniDir));

resDir = fullfile(iniDir, 'kwave_v1/simRSWF');
if ~exist("resDir","dir"); mkdir(resDir); end

simuNames = {'v1p6_homo'};
iSim = 1;

%% create the computational grid
Nx = 80;           % 120 number of grid points in the x (row) direction
Ny = 80;          % 120 number of grid points in the y (column) direction
Nz = 80;          % 120 page direction 
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
dz = 0.1e-3; 
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

%% Mask for inclusion
% OBS: In PAINT you can create a mask with desired inclusions shape
% save it as monochromatic file (.bmp) 
% mask Y.bmp is the name of the mask (Y shape)


% 
% maskInc=double(~imread("mask Y.bmp","bmp"));
% [x, y] = meshgrid(1:size(maskInc,2), 1:size(maskInc,1));
% 
% % Inclusion is repeated in the center of the cube (**centered in the X axis**)
% [xq, yq] = meshgrid(linspace(1, size(maskInc, 2), Ny), ...% Ny because 2nd dim
%                     linspace(1, size(maskInc, 1), Nz));   % Nz because 3rd dim
% maskInc_new = interp2(x, y, maskInc, xq, yq, "nearest");


%% MASK SIMPLE

centerLat = Ny / 2;  % lateral 
centerDep = Nx /2;   % depth
radius = Nx / 10;

% OPTION 1 
[Xgrid, Ygrid] = meshgrid(1:Ny, 1:Nx);
circ_mask = (Xgrid - centerLat).^2 + (Ygrid - centerDep).^2 <= radius^2;   

% OPTION 2 makeDisk
circ_mask2 = makeDisc(Nx, Ny, centerDep, centerLat, radius); 
circ_mask2 = logical(circ_mask2);

% OPTION 3 positions
ax_dep = kgrid.x(:,1,1); % in SI [m]
ax_lat = kgrid.y(1,:,1); % in SI [m]

centerLatm = 0e-3; % in SI [m]
centerDepm = 0e-3; % in SI [m]
radiusm = 1e-3;    % in SI [m]

circ_mask3m = (ax_lat-centerLatm).^2 + (ax_dep-centerDepm).^2 <= radiusm^2;
figure, imagesc(ax_lat*1e3, ax_dep*1e3, circ_mask3m)

% version homo
circ_mask3m = ones(size(circ_mask3m));

figure, 
subplot(121), imagesc(circ_mask)
subplot(122), imagesc(circ_mask2)
%
maskInc_new = circ_mask3m;
%% define the properties of the layers of the propagation medium

medium.sound_speed_shear       = zeros(Nx, Ny, Nz);    % [m/s]
medium.sound_speed_compression = zeros(Nx, Ny, Nz);   % [m/s]
medium.density                 = zeros(Nx, Ny, Nz);   % [kg/m^3]

% define the properties of the surrounding
medium.density                 (:, :, :) = 1000;     % [kg/m^3]
medium.sound_speed_shear       (:, :, :) = 1;      % [m/s]
medium.sound_speed_compression (:, :, :) = 1;   % [m/s] % fake compression wave speed in order to avoid divergence

thick=10; % (**centered in the X axis**) +/- thick = 2*thick+1 frames where the inclusion exist
medium.sound_speed_shear       ((0.5*Nx-thick):(0.5*Nx+thick),:,:) = 1+1*permute(repmat(maskInc_new,[1,1,2*thick+1]),[3,2,1]);      % [m/s]
medium.sound_speed_compression ((0.5*Nx-thick):(0.5*Nx+thick),:,:) = 1+1*permute(repmat(maskInc_new,[1,1,2*thick+1]),[3,2,1]);      % [m/s]
% fake compression wave speed in order to avoid divergence

% define the absorption properties
 medium.alpha_coeff_compression = 0.001; % [dB/(MHz^2 cm)]
 medium.alpha_coeff_shear       = 0.005; % [dB/(MHz^2 cm)]

% create the time array
%cfl   = 0.3;    % Courant-Friedrichs-Lewy number
%freq = 1000;
%t_end = 5/freq + kgrid.z_size/min(medium.sound_speed_shear(:));
% dTime = 1/(50*freq);
% kgrid.makeTime(max(medium.sound_speed_shear(:)), cfl, t_end);

dTime = 5e-6;
kgrid.setTime(3000,dTime)

% % VISUALIZATION €MZ
% figure,
% subplot(121), imagesc(medium.sound_speed_shear(:,:,1)), title('Speed shear')
% subplot(122), imagesc(medium.sound_speed_compression(:,:,1)),title('Speed compression')
%% 3D Plot and save simulated medium 

% Bottom half
xx = linspace(0,kgrid.x_size,Nx); % [m]
yy = linspace(0,kgrid.y_size,Ny); % [m]
zz = linspace(0,kgrid.z_size,Nz); % [m]
[~,~,Z] = meshgrid(xx,yy,zz);
mask3D_Cut=ones(size(medium.sound_speed_shear));
mask3D_Cut(Z<zz(Nz/2)) = 0;

fig = figure;
vol3d('CData',medium.sound_speed_shear,...
    'XData',1e2*[xx(1) xx(end)],...% [cm]
    'YData',1e2*[yy(1) yy(end)],...% [cm]
    'ZData',1e2*[zz(1) zz(end)],...% [cm]
    'Alpha',mask3D_Cut,...
    'texture','3D');
set(gca,'zdir','reverse');
xlabel('x (cm)'); ylabel('y (cm)'); zlabel('z (cm)');
axis tight; 
clim([0.5 2.5])
xlim(1e2*[xx(1) xx(end)]); ylim(1e2*[yy(1) yy(end)]); zlim(1e2*[zz(1) zz(end)]);

title('SWS');
view([-45,30]); grid on;
colormap('jet')
colorbar;
% saveas(fig,'SWS simulated medium - Bottom half.jpg');

saveas(fig, fullfile(resDir,"SWS_simMedium_BottomHalf_"+ simuNames{iSim}+".png"));

% Inclusion
xx = linspace(0,kgrid.x_size,Nx); %% m
yy = linspace(0,kgrid.y_size,Ny);
zz = linspace(0,kgrid.z_size,Nz);
[X,Y,Z] = meshgrid(xx,yy,zz);
mask3D_Cut=ones(size(medium.sound_speed_shear));
mask3D_Cut(1:(0.5*Nx-thick),:,:) = 0;
mask3D_Cut((0.5*Nx+thick):end,:,:) = 0;
                                
fig=figure;
Bmode_Processed = vol3d('CData',medium.sound_speed_shear,...
    'XData',1e2*[xx(1) xx(end)],...% [cm]
    'YData',1e2*[yy(1) yy(end)],...% [cm]
    'ZData',1e2*[zz(1) zz(end)],...% [cm]
    'Alpha',mask3D_Cut,...
    'texture','3D');
set(gca,'zdir','reverse');
xlabel('x (cm)'); ylabel('y (cm)'); zlabel('z (cm)');
axis tight; 
clim([0.5 2.5])
xlim(1e2*[xx(1) xx(end)]); ylim(1e2*[yy(1) yy(end)]); zlim(1e2*[zz(1) zz(end)]);

title('SWS');
view([-45,30]); grid on;
colormap('jet')
colorbar;
% saveas(fig,'SWS simulated medium - Inclusion.jpg');


saveas(fig, fullfile(resDir,"SWS_simMedium_Inc_"+ simuNames{iSim}+".png"));

%%  source
source.u_mask = zeros(Nx,Ny,Nz);
numSources = 6;

% Get the surface indices
surface_indices = [];
% Front and back faces
[ix, iy] = ndgrid(1:Nx, 1:Ny);
surface_indices = [surface_indices; [ix(:), iy(:), ones(Nx*Ny,1)]; ...       % Front face
                              [ix(:), iy(:), Nz*ones(Nx*Ny,1)]];              % Back face

% Left and right faces
[iz, iy] = ndgrid(1:Nz, 1:Ny);
surface_indices = [surface_indices; [ones(Nz*Ny,1), iy(:), iz(:)]; ...       % Left face
                              [Nx*ones(Nz*Ny,1), iy(:), iz(:)]];              % Right face

% Top and bottom faces
[ix, iz] = ndgrid(1:Nx, 1:Nz);
surface_indices = [surface_indices; [ix(:), ones(Nx*Nz,1), iz(:)]; ...       % Top face
                              [ix(:), Ny*ones(Nx*Nz,1), iz(:)]];              % Bottom face

% Randomly select points from the surface
num_surface_points = size(surface_indices, 1);
random_indices = randperm(num_surface_points, numSources);
selected_points = surface_indices(random_indices, :);

% Assign 1 to the selected random surface points
for i = 1:numSources
    source.u_mask(selected_points(i,1), selected_points(i,2), selected_points(i,3)) = 1;
end

source_freq = 500;
source_mag = 1;
source.uz = source_mag*sin(2*pi*source_freq*kgrid.t_array);

%% Define a sensor
sensor.mask = mask3D_Cut;
sensor.record = {'u','u_split_field'};

%%
% define input arguments
% input_args = {'PlotScale', [-0.75, 0.75, -0.15, 0.15], 'PlotPML', false,...
%      'DisplayMask', display_mask, 'DataCast', 'single'};
input_args = {'DataCast',DATA_CAST};

sensor_data = pstdElastic3D(kgrid, medium, source, sensor, input_args{:});

%% save Attempt 1
myfileName = ['SensorData_RSWF_',simuNames{iSim},'.mat'];

% save all (DEPRECATED
% save(fullfile(resDir, myfileName), 'sensor_data', 'kgrid', '-v7.3');


%% save Attempt 2
myfileName = ['RSWF_',simuNames{iSim},'.mat'];

uz_s = sensor_data.uz_split_s;

% modify according to the size of the sensor (Nx,..,Nz,frames);
uz_s = reshape(uz_s,2*thick-1,Ny,Nz,[]);

% save('Planar4Dfull.mat', 'Planar4Dfull50','-v7.3')

uz_s_array = gather(uz_s);
save(fullfile(resDir, myfileName), 'uz_s_array', 'kgrid', '-v7.3');

