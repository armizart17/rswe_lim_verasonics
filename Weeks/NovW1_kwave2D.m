% SIMPLE SCRIPT FOR 2D SIMULATION

clc; clear; close all;

addpath(genpath('/opt/MATLAB Add-Ons/')); 
addpath(genpath(pwd));

disp(pwd);

DATA_CAST = 'gpuArray-single';  % GPU
% DATA_CAST = 'single';           % CPU

iniDir = '/mnt/nfs/emiranda/proj/rswe/'; % for cluster

addpath(genpath(iniDir));

resDir = fullfile(iniDir, 'kwave_v1/simRSWF_2D');
if ~exist("resDir","dir"); mkdir(resDir); end

simuNames = {'v1'};
iSim = 1;

%% create the computational grid for 2D
Nx = 80;           % number of grid points in the x (row) direction
Ny = 80;           % number of grid points in the y (column) direction
dx = 0.1e-3;       % grid point spacing in the x direction [m]
dy = 0.1e-3;       % grid point spacing in the y direction [m]

% Define the 2D computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% Mask for inclusion (2D mask)
centerLat = Ny / 2;  % lateral 
centerDep = Nx /2;   % depth
radius = Nx / 10;

% Create a 2D circular mask
[Xgrid, Ygrid] = meshgrid(1:Ny, 1:Nx);
circ_mask = (Xgrid - centerLat).^2 + (Ygrid - centerDep).^2 <= radius^2;   

figure, imagesc(circ_mask), title('Circular Mask for Inclusion')

%% define the properties of the layers of the propagation medium

medium.sound_speed_shear       = zeros(Nx, Ny);    % [m/s]
medium.sound_speed_compression = zeros(Nx, Ny);    % [m/s]
medium.density                 = zeros(Nx, Ny);    % [kg/m^3]

% define the properties of the surrounding
medium.density(:,:)                 = 1000;       % [kg/m^3]
medium.sound_speed_shear(:,:)       = 1;          % [m/s]
medium.sound_speed_compression(:,:) = 1;          % [m/s] % fake compression wave speed to avoid divergence

% Define the inclusion
medium.sound_speed_shear(circ_mask)       = 2;    % [m/s]
medium.sound_speed_compression(circ_mask) = 2;    % [m/s]

% Define the absorption properties
medium.alpha_coeff_compression = 0.001;  % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 0.005;  % [dB/(MHz^2 cm)]

% Create the time array
dTime = 5e-6;
kgrid.setTime(3000, dTime);

% Visualization of the medium
figure,
subplot(121), imagesc(medium.sound_speed_shear), title('Shear Speed (2D)')
subplot(122), imagesc(medium.sound_speed_compression), title('Compression Speed (2D)')

%% Define a 2D source
source.u_mask = zeros(Nx, Ny);
numSources = 3;

% Randomly assign sources within the grid (front face only for simplicity)
source_positions = randi([1 Nx*Ny], 1, numSources); % Random source positions
source.u_mask(source_positions) = 1;

source_freq = 500; % [Hz]
source_mag = 1;
source.ux = source_mag * sin(2 * pi * source_freq * kgrid.t_array);  % X-component of displacement

%% Define a 2D sensor
sensor.mask = circ_mask;
sensor.record = {'u', 'u_split_field'};

%% Define input arguments
input_args = {'DataCast', DATA_CAST};

% Run the simulation for 2D
sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

%% Save the sensor data
myfileName = ['SensorData_RSWF_2D_', simuNames{iSim}, '.mat'];
save(fullfile(resDir, myfileName), 'sensor_data', 'kgrid', '-v7.3');

