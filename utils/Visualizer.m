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
    end
end
