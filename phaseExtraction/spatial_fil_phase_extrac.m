function [Frames0,Frames1] = spatial_fil_phase_extrac(u_new, peak, L)
% function [Frames0,Frames1] = spatial_fil_phase_extrac(u_new, peak, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After temporal fintering, a 2-D FIR filter using 2-D window method is performed 
% using a % spatial bandpass filter with cutoffs cs_min & cs_max. 
% Then 4 complex matrixes are obtained: 
% Inputs:      
%         f           : vibration frequency
%         peak        : where the maximun spectrum point is  located        
%         u_new       : temporal filtered particle velocity data  
%         cs_min(max) : cutoffs frequencies for bandpass spatial filter
% Outputs: 
%         Frames0       : magnitude and phase data before spatial filtering
%         Frames1       : only phase data before spatial filtering
%         Frames2       : magnitude and phase data after spatial filtering
%         Frames3       : only phase data after filtering
%%%%% Results using Frames0 & Frames3 are presented in Ormachea et al (2018)
% Author: Edited by EMZ from LIM repository
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sigma = 50; 
% Fs = 1/dinf.dx;                         % Sampling spatial frequency

[~,Ps,Frames0] = amp_phs(u_new,peak,L); % magnitude and phase data before spatial filtering

% Vw1 = fftshift(fft2(ifftshift(Frames0))); % spactial frequency domain for filtering

% Spatial frequencies cutoffs estimation base on  the relationship  k=2pi/c 

% [k1,k2] = freqspace(size(Vw1),'meshgrid');
% Hd = ones(size(Vw1)); 
% r = sqrt(k1.^2 + k2.^2);
% kl = (2*pi/cs_max)*f*2/(2*pi*Fs); kh = (2*pi/cs_min)*f*2/(2*pi*Fs);
% Hd((r<kl)|(r>kh)) = 0;
% win = fspecial('gaussian',size(Vw1),sigma); 
% win = win ./ max(win(:));  % Make the maximum window value be 1.
% h = fwind2(Hd,win);        % Using the 2-D window, design the filter that best produces the desired frequency response
% mask2 = abs(real(fftshift(fft2(ifftshift(h)))));

 
% new_Vw1 = Vw1.*mask2;                          % mask filtering in spatial frequency domain 
% un2 = zeros(size(u_new));
% for ii=1:size(u_new,3)
%     un1 = u_new(:,:,ii);
%     un2(:,:,ii) = filter2(h,un1);
% end

% [~,Ps1,Frames2] = amp_phs(un2,peak,L); % magnitude and phase data after spatial filtering

% Frames2 = fftshift(ifft2(ifftshift(new_Vw1))); % magnitude and phase data after spatial filtering
% Ps1 = angle(Frames2);

Frames1 = exp(1i*Ps);  %  only phase data before spatial filtering
% Frames3 = exp(1i*Ps1);    %  only phase data after filtering

% Frames0 = magnitude*exp(jPhase) : magnitude & phase from temporal data
% Frames1 = 1*exp(jPhase)
% Frames2 = magnitude2*exp(jPhase2): magnitude2 & phase2 after spatial filt
% Frames3 = 1*exp(jPhase2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Uncomment if you want to plot Magnitude 2D complex matrix, filter profiles, 2D-mask, magnitude and
%%%%%%%% phase before and after spatial filtering
% 
% res_x=dinf.dx;
% res_y=dinf.dz;
% [m,n,~]=size(Frames0);
% x_vec=0:res_x:res_x*(n-1);
% y_vec=0:res_y:res_y*(m-1);
% y_vec = y_vec+dinf.u_rows(1)*dinf.dz;
% 

end