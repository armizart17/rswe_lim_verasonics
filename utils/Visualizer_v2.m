classdef Visualizer_v2
    % Class to handle visualization of ultrasound data with imOverlayInterp
    
    properties
        x;           % Lateral coordinates
        z;           % Axial coordinates
        BmodeFull;   % B-mode image data
        x_roi;       % ROI lateral coordinates
        z_roi;       % ROI axial coordinates
        img_big;     % Image to overlay
        caxis_bmode; % Color axis for B-mode
        caxis_img;   % Color axis for overlay Image
        title_name;  % Title for the plot
        fact_transparency; % Transparency factor
        units;       % Units for plotting ('cm' or 'mm')
    end
    
    methods
        function obj = Visualizer_v2(defaultStruct)
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
            % Calculate the mean of img_big
            mean_sld = mean(obj.img_big(:), 'omitnan');
        end
        
        function sig_sld = std(obj)
            % Calculate the standard deviation of img_big
            sig_sld = std(obj.img_big(:), 'omitnan');
        end
        
        function visualize(obj)
            % Method to visualize the data using imOverlayInterp
            
            % Title format
            pmSymbol = char(177); % ±
            title_str = sprintf('%s: %.2f %c %.2f', obj.title_name, obj.mean(), pmSymbol, obj.std());

            % Scale factors
            if strcmp(obj.units, 'mm')
                scaleFactor = 1e3; % Convert to mm
                unit = '[mm]';
            else
                scaleFactor = 1e2; % Convert to cm
                unit = '[cm]';
            end

            % Create ROI mask (Assume all values inside `x_roi` and `z_roi` are ROI)
            [X,Z] = meshgrid(obj.x, obj.z);
            ROI = (X >= min(obj.x_roi) & X <= max(obj.x_roi)) & ...
                  (Z >= min(obj.z_roi) & Z <= max(obj.z_roi));

            % Call the overlay function
            figure;
            set(gcf, 'Position', [10, 10, 600, 700]);

            [~,~,hColor] = imOverlayInterp(...
                obj.BmodeFull, obj.img_big, ...
                obj.caxis_bmode, obj.caxis_img, ...
                obj.fact_transparency, ...
                scaleFactor * obj.x, scaleFactor * obj.z, ...
                ROI, scaleFactor * obj.x_roi, scaleFactor * obj.z_roi);

            % Formatting
            xlabel(['\bfLateral ', unit], 'fontsize', 14);
            ylabel(['\bfDepth ', unit], 'fontsize', 14);
            title(title_str, 'Interpreter', 'none', 'FontWeight', 'bold', 'fontsize', 14);
            set(gca, 'fontsize', 14);

            % Adjust colorbar
            ylabel(hColor, 'm/s', 'fontsize', 14);
            colormap jet;
        end
    end
end
