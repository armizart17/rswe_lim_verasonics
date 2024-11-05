function create_pv_movie(u,dinf,framerate,freq,filtered)
% function create_pv_movie(u,dinf,framerate,freq,filtered)
% Create a video writer object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Video maker to observe particle velocity (pattern of waves)
% 
% Inputs:  u         - Matrix containing the particle velocity after 
%                      Loupas' estimation (indicate if filtered in last 
%                      arg of this function)
%          dinf      - Structure that contains at least dx, dz and Tprf
%          framerate - framerate of the desired video
%          freq      - Vibrational frequency of harmonic waves
%          filtered  - 0 = unfiltered pv
%                      1 = filtered pv
% Outputs: .mp4 file saved on the current directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if filtered == 0
    aux='unfilt';
elseif filtered == 1
    aux='filt';
else
    error('Filtered should be 0 or 1');
end
    outputVideo = VideoWriter([aux '_part_veloc_' num2str(freq) 'Hz.mp4'], 'MPEG-4');
    outputVideo.FrameRate = framerate;  % Set the frame rate
    
    % Open the video writer object to start writing frames
    open(outputVideo);
    
    % Vector for axis and Number of frames
    % xvec=dinf.dx*(0:size(u,2)-1)*1e2;
    % zvec=dinf.dz*(0:size(u,2)-1)*1e2;
    
    % MODIFICATION EMZ
    xvec = dinf.dx*((-(size(u,2)-1)/2) : (size(u,2)-1)/2) * 1e2; % [cm]
    zvec = dinf.dz*(0:size(u,1)-1) * 1e2;                        % [cm] 
    nFrames = size(u,3);
    
    % Create the figure
    hfig=figure;
    himag= imagesc(xvec,zvec,u(:,:,1));
    %caxis(
    xlim([min(xvec) max(xvec)]);
    ylim([min(zvec) max(zvec)]);
    xlabel('Lateral (cm)','FontName', 'Cambria', 'FontSize', 18);
    ylabel('Depth (cm)','FontName', 'Cambria', 'FontSize', 18);
    
    
    % Loop through to create figures and write them to the video
    for kk = 2:nFrames

       set(himag,'CData',u(:,:,kk));  
       title(sprintf('Frame %d, time = %0.2f ms', kk, dinf.prf*1e3*kk),'FontName', 'Cambria', 'FontSize', 18);
       drawnow;
        
        % Capture the current figure as a frame
        frame = getframe(hfig);  
        
        % Write the frame to the video
        writeVideo(outputVideo, frame);

    end
    
    % Close the video writer object to finish the video file
    close(outputVideo);
    close(hfig);
    
    % Display a message when done
    disp('Video creation complete!');
end