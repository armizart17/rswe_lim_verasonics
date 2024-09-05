function [fu1, Wave_z, Frames1] = u_filt(u_new, f, f_band, dinf, cs_min, cs_max)
% function [fu1, Wave_z, Frames1] = u_filt(u_new, f, f_band, dinf, cs_min, cs_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal filtering process of particle velocity signals: a median filter
% to reduce noise and peppper noise and the a bandpass Hamming FIR filter are applied
% Inputs:  
%           u_new   : input particle velocity
%           f       : Vibration freqeuncy
%           f_band  : Frequency range for the bandpass cuttoffs are +-
%                    2*f_band
%           u_new   : particle velocity data with 10 periods
% Outputs: 
%           fu1     : wave  
%           wave_z  : magnitude and phase data after spatial filtering   
%           Frames1 : only phase data after spatial filtering
% Author: Edited by EMZ from LIM repository
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Fs = dinf.PRFe;  % Sampling frequency
    %$  A 50th-order FIR bandpass filter with bandpass [fv-10 dv+10]/Fs cutoffs. 
    ord = 200; % filter order
    b = fir1(ord,[2*(f-f_band)/Fs 2*(f+f_band)/Fs],tukeywin(ord+1,1)); 
    fu = zeros(size(u_new));
    for ii=1:size(fu,1)
        for jj=1:size(fu,2)
            s1 = squeeze(u_new(ii,jj,:));
            fu(ii,jj,:) = filter(b,1,s1);
        end
    end
    
    
    %%
    
    Fs = dinf.PRFe;                        % Sampling frequency
    L = 1e4;                               % Length of signal
    df = Fs/L;                              % frequency resolution
    f1 = (-round(L/2):1:round(L/2)-1)*df;   % frequency axis
    [ ~, ix ] = min( abs( f1-f ));  % find index for closest frequency to the vibration frequency
    peak = ix;
    [Frames0,~] = spatial_fil_phase_extrac(fu,peak,L);
    
    %%
    sigma = 300; 
    Fs1 = 1/dinf.dz;    % Sampling spatial frequency
    Fs2 = 1/dinf.dx;
    
    % Vw1 = squeeze(fu(:,:,5));
    % Spatial frequencies cutoffs estimation base on  the relationship  k=2pi/c 
    [k1,k2] = freqspace(size(Frames0),'meshgrid');
    k1 = k1*(2*pi*Fs1);
    k2 = k2*(2*pi*Fs2);
    
    Hd = ones(size(Frames0)); 
    r = sqrt(k1.^2 + k2.^2);
    kl = (2*pi*f/cs_max)*2; kh = (2*pi*f/cs_min)*2;
    
    Hd((r<kl)|(r>kh)) = 0;
    win = fspecial('gaussian',size(Frames0),sigma); 
    win = win ./ max(win(:));  % Make the maximum window value be 1.
    h = fwind2(Hd,win);        % Using the 2-D window, design the filter that best produces the desired frequency response
    % mask2 = abs(real(fftshift(fft2(ifftshift(h)))));
    
     
    % mask filtering in spatial frequency domain 
    % for ii=1:size(fu,3)
    %     un1 = fu(:,:,ii);
    %     un1 = medfilt2(abs(Frames0),[11 7]);
        Wave_z = filter2(h,Frames0);
    % end
    
    omega = 2*pi*f;
    resT = 1/Fs;  
    Tmax = 60*1e-3;  
    t = 0:resT:Tmax;
    fu1 = zeros([size(Wave_z),length(t)]);
    for kk = 1:length(t)
        fu1(:,:,kk) = abs(Wave_z).*cos(angle(Wave_z)+omega*t(kk));
    end
    
    Frames1 = exp(1i*angle(Wave_z));
end

%%
% figure, 
%     for kk = 1:length(t)
%         fu1(:,:,kk) = abs(Wave_z).*cos(angle(Wave_z)+omega*t(kk));
%         imagesc(real(fu1(:,:,kk)))
%         pause(0.01);
%     end
% 
%     %%
% figure, 
% imagesc(angle(Frames1));
% 
% figure, 
% imagesc(angle(Wave_z));
% 
% figure, 
% imagesc(angle(Frames1) - angle(Wave_z))
% 
% isequal(angle(Frames1), angle(Wave_z))

%%
%% EMZ
    % [h, w] = freqz(b, 1, [], Fs);
    % figure;
    % plot(w, mag2db( abs(h) )), hold on, grid minor,
    % xline(f-f_band, 'k--');
    % xline(f+f_band, 'k--');
    % xlabel('Frequency (Hz)');
    % ylabel('Magnitude');
    % title('Magnitude Response in Hz');
    % 
    % 
    % ii = ceil(size(u_new, 1) / 2);
    % jj = ceil(size(u_new, 2) / 2);
    % 
    % figure(17), 
    % set(gcf, 'units','normalized', "Position", [0 0 720 980])
    % for ii = 1:size(u_new, 1)
    %     for jj = 1:size(u_new, 2)
    % 
    %         signal_og = squeeze( squeeze (u_new(ii, jj, :) ) );
    %         signal_fi = squeeze( squeeze (fu(ii, jj, :) ) );
    % 
    %         [fHz, ~, Y_og] = fft_info(signal_og, Fs);
    %         [fHz, ~, Y_fi] = fft_info(signal_fi, Fs);
    % 
    % 
    %         subplot(2,1,1), plot(fHz, mag2db( abs(Y_og) ),'k', 'DisplayName', 'Og'), grid on;
    %         legend;
    %         xlim([-1000 1000])
    %         subplot(2,1,2), plot(fHz, mag2db( abs(Y_fi) ),'r' , 'DisplayName', 'Fi'), grid on;
    %         xlim([-1000 1000]);
    %         xline(f-f_band, 'k--');
    %         xline(f+f_band, 'k--');
    %         % 
    %         % subplot(2,1,1), plot(fHz, mag2db( abs(Y_og) ),'k', 'DisplayName', 'Og'), grid on;
    %         % legend;
    %         % subplot(2,1,2), plot(fHz, mag2db( abs(Y_fi) ),'r' , 'DisplayName', 'Fi'), grid on;
    %         % xlabel('Freq [Hz]'), ylabel('Magn')
    %         legend;
    %         pause(0.01)
    % 
    %     end
    % 
    % end
    % 
