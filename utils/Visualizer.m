%% CLASS FOR VISUALIZATION  %%
% Description: 
% A generic class for Image Visualization (Superposition with Bmode)
% Version 2.0 (Designed for Color Imaging) 
% Author: EMZ (based on LIM codes)

classdef Visualizer
    % Class to handle visualization of ultrasound data
    
    properties
        x;           % Lateral coordinates
        z;           % Axial coordinates
        BmodeFull;   % B-mode image data
        x_roi;       % ROI lateral coordinates
        z_roi;       % ROI axial coordinates
        img_big;     % Imagen big usually after interp2 or resize
        caxis_bmode; % Color axis for B-mode
        caxis_img;   % Color axis for Image
        title_name;  % Title for the plot
        fact_transparency; % Transparency factor
        units;       % Units for plotting ('cm' or 'mm')
    end
    
    methods
        function obj = Visualizer(defaultStruct)
            % Constructor to initialize the object with data
            obj.x = defaultStruct.x;
            obj.z = defaultStruct.z;
            obj.BmodeFull = defaultStruct.BmodeFull;
            obj.caxis_bmode = defaultStruct.caxis_bmode;
            obj.caxis_img = defaultStruct.caxis_img;
            obj.fact_transparency = defaultStruct.fact_transparency;

        end

        function obj = setUnits(obj, units)
            % Setter for units with validation
            if strcmp(units, 'cm') || strcmp(units, 'mm')
                obj.units = units;
            else
                error('Units must be either "cm" or "mm".');
            end
        end

        function obj = setROI(obj, x_roi, z_roi, img_big)
            obj.x_roi = x_roi;
            obj.z_roi = z_roi;
            obj.img_big = img_big;
        end

        function obj = setTitle(obj, title_name)
            obj.title_name = title_name;
        end

        function mean_sld = mean(obj)
            % Method to calculate the mean of ACS_big
            mean_sld = mean(obj.img_big(:));
        end
        
        function sig_sld = std(obj)
            % Method to calculate the standard deviation of ACS_big
            sig_sld = std(obj.img_big(:));
        end
        
        function visualize(obj)
            % Method to visualize the data
                        
            % Create title string
            pmSymbol = char(177); %+/-
            title_str = sprintf('%s: %.2f %c %.2f', obj.title_name, obj.mean(), pmSymbol, obj.std());
            
            % Determine scale factor and unit label
            if strcmp(obj.units, 'mm')
                scaleFactor = 1e3; % Convert from cm to mm
                unit = '[mm]';
            else
                scaleFactor = 1e2;  % Default to cm
                unit = '[cm]';
            end
            
            % Plot B-mode image
            figure;
            set(gcf, 'Position', [10 10 600 700])
            subimage(scaleFactor* obj.x, scaleFactor* obj.z, ...
                     256 * mat2gray(obj.BmodeFull, obj.caxis_bmode), gray(256));
            hold on;
            
            % Plot ACS_sld_big with transparency
            h = subimage(scaleFactor* obj.x_roi, scaleFactor*  obj.z_roi, ...
                         256 * mat2gray(obj.img_big, obj.caxis_img), jet(256));
            set(h, 'AlphaData', 1 - obj.fact_transparency);
            axis('tight')
            
            % Set labels and title
            xlabel(['\bfLateral ', unit], 'fontsize', 14);
            ylabel(['\bfDepth ', unit], 'fontsize', 14);
            title(title_str, 'Interpreter', 'none', 'FontWeight', 'bold', 'fontsize', 14);
            set(gca, 'fontsize', 14);
            
            % Add colorbar
            h2 = colorbar;
            ylabel(h2, 'm/s','fontsize', 14);
            caxis(obj.caxis_img);
            colormap(jet);
        end

        function [rect_scaled, bin_mask] = selectMetricRec(obj)
            % Method to select ROI and calculate mean and std of selected region
            
            % Determine scale factor based on units
            if strcmp(obj.units, 'mm')
                scaleFactor = 1e3; % Convert to mm
                unitLabel = '[mm]';
            else
                scaleFactor = 1e2; % Convert to cm
                unitLabel = '[cm]';
            end

            % Display the image
            lw = 2;
            figure, 
            imagesc(scaleFactor * obj.x, scaleFactor * obj.z, obj.img_big, obj.caxis_img);
            colorbar, colormap(jet(256))
            xlabel(['\bfLateral ', unitLabel]);
            ylabel(['\bfDepth ', unitLabel]);
            title('Image map');
            
            % Select the ROI using getrect
            confirmation = '';
            while ~strcmp(confirmation,'Yes')
                rect = getrect; % Allow the user to select the rectangular ROI
                confirmation = questdlg('Sure?');
                if strcmp(confirmation,'Cancel')
                    disp(rect);
                    break
                end
            end
            close,

            % Convert the selected rectangle back to meters (if in mm or cm)
            rect_scaled = rect / scaleFactor; % To [m]

            x_inf = rect_scaled(1); 
            x_sup = rect_scaled(1) + rect_scaled(3);
            z_inf = rect_scaled(2); 
            z_sup = rect_scaled(2) + rect_scaled(4);

            % Limit data to the selected region
            ind_x = (x_inf <= obj.x) & (obj.x <= x_sup);
            ind_z = (z_inf <= obj.z) & (obj.z <= z_sup);

            obj.x_roi = obj.x(ind_x);
            obj.z_roi = obj.z(ind_z);

                % Create a binary mask based on the selected ROI
            [X, Z] = meshgrid(obj.x, obj.z); % Create a grid of x and z coordinates
            bin_mask = (X >= x_inf & X <= x_sup) & (Z >= z_inf & Z <= z_sup);
            % image_region = obj.img_big(mask); % Extract the 2D region

            image_region = obj.img_big(ind_z, ind_x); % Extract the 2D region
            image_region_1D = image_region(:); % Flatten the 2D matrix to 1D

            % Calculate mean and standard deviation
            mean_reg = mean(image_region_1D);
            std_reg = std(image_region_1D);

            % Display the selected region with a rectangle overlay
            pm_char = char(177); % The ± symbol
            title_str = sprintf('Reg. Metrics: %.2f %c %.2f', mean_reg, pm_char, std_reg);

            figure, 
            imagesc(scaleFactor * obj.x, scaleFactor * obj.z, obj.img_big, obj.caxis_img), hold on;
            rectangle('Position', rect, 'EdgeColor', 'k', 'LineWidth', lw, 'LineStyle', '-.');
            xlabel(['\bfLateral ', unitLabel]);
            ylabel(['\bfDepth ', unitLabel]);
            title(title_str);
            
            colorbar, colormap(jet(256))

        end
    end
end
