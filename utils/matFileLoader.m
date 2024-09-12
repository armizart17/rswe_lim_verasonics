function typName = matFileLoader
    % Create a figure for the GUI
    f = figure('Name', 'MAT File Loader', ...
               'Position', [300, 300, 400, 150], ...
               'MenuBar', 'none', ...
               'NumberTitle', 'off', ...
               'WindowStyle', 'modal'); % Modal to block script execution
    
    % Create a button for loading a .mat file
    uicontrol('Style', 'pushbutton', ...
              'String', 'Load .mat file', ...
              'Position', [150, 70, 120, 60], ...
              'FontSize', 12, ...
              'HorizontalAlignment', 'center', ...
              'Callback', @loadMatFile);

    % Create a text area to display the file name
    fileText = uicontrol('Style', 'text', ...
                         'Position', [50, 20, 300, 30], ...
                         'String', 'No file loaded', ...
                         'FontSize', 11, ...
                         'HorizontalAlignment', 'center');

    % Predefine the typName output as empty
    typName = '';

    % Inner function to load the .mat file
    function loadMatFile(~, ~)
        % Open file dialog to select a .mat file
        [fileName, pathName] = uigetfile('*.mat', 'Select a .mat file');
        
        % If a file is selected, proceed to load the .mat file
        if fileName
            fullname = fullfile(pathName, fileName); % Full file path
            
            % Load the variables from the .mat file into the base workspace
            data = load(fullname);
            
            % Get the typName (filename without extension)
            typName = fileName(1:end-4);  % Remove '.mat'
            
            % Import all variables from the .mat file into the base workspace
            varNames = fieldnames(data);
            for i = 1:numel(varNames)
                assignin('base', varNames{i}, data.(varNames{i}));
            end
            
            % Update the text area to display the loaded file name
            fileText.String = ['Loaded: ', fileName];
            
            % Close the GUI after successful loading
            uiresume(f); % Resume UI execution to return control
            close(f);
        else
            % If no file is selected, update the text area to notify the user
            fileText.String = 'No file loaded';
        end
    end
    
    % Wait for the GUI to load a file before returning
    uiwait(f);
end
