%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix creation with complex data oscillating at the vibration frequency
% in the temporal domain. No need for temporal filtering since the
% vibration frequency is already known.
% 
% Inputs:  u           - Particle velocity signal
%          f_v         - Vibration freqeuncy
%          dinf        - Structure that contains the PRF
%
% Outputs: u_complex  - Complex 2-D matrix with desired data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_complex = temporal_extraction(u,f_v,Fs)
    %freq_range = Fs*(0:ceil(size(u,3)/2))/ceil(size(u,3));
    n = 2^nextpow2(size(u,3)*1  );
    freq_range = linspace(-round(Fs/2),round(Fs/2),n);
    index = find_position(f_v,freq_range);
    size_matrix = [size(u,1),size(u,2)];
    u_complex = zeros(size_matrix);
    figure, 
    for ii=1:size(u,1)
        for jj=1:size(u,2)
            u_time = squeeze(u(ii,jj,:))';
            %plot(t_vector,u_TempFilt)
            u_freq = fftshift(fft(u_time,n));
            if ii == floor(size(u,1)/2) && jj == floor(size(u,2)/2)
                %scale_x=[0 260  800 1200 2000];
                
                plot(freq_range,20*log10( abs(u_freq)./abs(max(u_freq)) ),'k','Linewidth',1); 
                xlabel('Frequency [Hz]'), ylabel('Normalized Magnitude')
                title('Frequency of raw data'), grid on
                % set(gca,'Fontname','Times')
                %set(gca,'XTick',scale_x,'XtickLabel',scale_x);
            end
            if ~isempty(index)
                u_complex(ii,jj) = u_freq(index); % extracts the phase 
            else
                u_complex(ii,jj) = 0;
            end
        end
    end
end