%% Script to plot l2-norm of gradient of each kernel for visualization
% EMZ, Oct 2024
%%
% Display the image using imagesc
figure, 
imagesc(angle(pv_field_pad));
colorbar;  % Optional: Add a colorbar
hold on;   % Keep the plot to overlay the rectangle

% Define the position and size of the rectangle
xPos = search_x(jj);  % X-coordinate (horizontal position)
yPos = search_z(ii);  % Y-coordinate (vertical position)
width = w_kernel(2);  % Width of the rectangle (15 pixels)
height = w_kernel(1);  % Height of the rectangle (15 pixels)

% Draw the rectangle
rectangle('Position', [xPos, yPos, width, height], 'EdgeColor', 'r', 'LineWidth', 2);

% Optional: Add a title to the figure
title('Phase with 15x15 kernel');

hold off;  % Release the hold on the plot
%%
% Parallel loop over rotations using the given number of angles
clear val_abs_3D
num_angles = 6;
    for idx = 1:num_angles

        clear diff_area_z diff_area_x
        % Calculate the actual rotation angle
        theta = mod((idx - 1) * delta_theta, 180); % Ensuring the angle wraps around at 360

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
%%


% Number of channels
numChannels = size(val_abs_3D, 3);

% Example list of angles


% Create subplots for each channel
figure;

nRows = 2;
nCols = ceil(numChannels/nRows);

tiledlayout(nRows, nCols)

for i = 1:numChannels
    nexttile
    myAngle = mod((i - 1) * delta_theta, 180);
    % subplot(4, 3, i);  % Create a 2x3 grid of subplots (since you have 6 channels)
    imagesc(val_abs_3D(:,:,i), [0 3000]);  % Plot the 2D slice for each channel
    colorbar;  % Optional: Add a colorbar to each subplot
    hold on;
    contour(circ_mask, [1 1], 'r-', 'LineWidth', 2); % Plot the contour of the binary mask

    title(['\theta = ' num2str(myAngle) '°']);  % Add title with the corresponding angle
    axis square;  % Make each subplot square
end

% Optional: Add a main title for the entire figure
sgtitle('Grad-l2 eq K');
%%

% Number of channels
numChannels = size(val_abs_3D, 3);

% Create subplots for each channel
figure;
tiledlayout(nRows, nCols)
nRows = 2;
nCols = ceil(numChannels/nRows);

for i = 1:numChannels

    % Calculate angle
    myAngle = mod((i - 1) * delta_theta, 180);

    % Create subplot in a 4x3 grid
    % subplot(4, 3, i);  
    nexttile
    
    % Create X and Y grids for surf
    [X, Y] = meshgrid(1:size(val_abs_3D, 2), 1:size(val_abs_3D, 1));

    % Use surf instead of imagesc to plot 3D surface
    surf(X, Y, val_abs_3D(:,:,i));
    
    % Optional: Set z-axis limits and colormap scaling
    zlim([0 3000]);
    caxis([0 3000]);  % Adjust colormap scaling to match values
    colorbar;  % Add a colorbar to each subplot
    
    % Hold to overlay contour
    hold on;
    contour(X, Y, circ_mask, [1 1], 'r-', 'LineWidth', 2); % Plot the contour of the binary mask

    % Set the title with the angle
    title(['\theta = ' num2str(myAngle) '°']);  
    
    % Adjust the appearance
    axis square;
    shading interp;  % Smooth shading for the surface
end

% Optional: Add a main title for the entire figure
sgtitle('Grad-l2 eq K');