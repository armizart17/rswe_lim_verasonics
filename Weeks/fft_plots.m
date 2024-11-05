    %clear;close all; clc;
[u,dinf] = pv_cal(IQData1,dinf, dinf.num_angles); 
u = signal_period(freq, dinf.PRFe, u); 
%load('particle_displacement_500Hz_9percent.mat');
%%
prf=dopPRF; %Sampling frequency of the transductor in Hz
resol_x=dinf.dx; % resolution in m 
resol_z=dinf.dz; % resolution in m 

f_spacial_x=1/resol_x; % SpatialFrequency in 1/m
f_spacial_z=1/resol_z; % SpatialFrequency in 1/m
f_v=freq; % Vibration frequency of the external source
c_max=6; c_min=0.5; % maximun & minumin admissible SWS
k_min=(2*pi*f_v)/c_max; k_max=(2*pi*f_v)/c_min; % max & min admissible wavenumbers
freq_range=prf*(0:size(u,3)/2)/size(u,3); % Frequency scale depending on sampling frequency and size of the signal
% t_vector=(0:size(u,3)-1)*(1/prf); 

xx = xdim;
yy = ydim;

% Finding the index in which it's contained the closest value to the f_v
index_0=min(abs(f_v-freq_range)); 
index_1=find(freq_range==f_v+min(abs(f_v-freq_range)));
if index_1~=0
    index=index_1;
else
    index=find(freq_range==f_v-index_0);
end

%% Frequency filtering
%n_2=2^nextpow2(size(u,3));

coef=zeros(size(u,1),size(u,2));
for a=1:size(u,1)
    for b=1:size(u,2)
        u_TempFilt=transpose(squeeze(u(a,b,:)));
        %plot(t_vector,u_TempFilt)
        %u_freq=fft(u_TempFilt,n_2);
        u_freq=fft(u_TempFilt);
        if a==fix(size(u,1)/2) && b==fix(size(u,2)/2)
            scale_x=[0 100  400 1200 2000];
            fig1=figure;
            plot(freq_range,abs(u_freq(1:size(u,3)/2 +1))./abs(max(u_freq))); xlabel('Frequency [Hz]');ylabel('Normalized Magnitude');title('Frequency domain of raw data');
            set(gca,'XTick',scale_x,'XtickLabel',scale_x);
            % saveas(fig1,'FrequencyPeak.png');
        end
        if index~=0
            coef(a,b)=u_freq(index); % extracts the phase 
        else
            coef(a,b)=0;
        end
    end
end

%% PLOT 9 POINTS 

coef=zeros(size(u,1),size(u,2));
count = 0; % Counter for subplots
fig = figure; % Create a new figure for subplots

for a = 1:size(u,1)
    for b = 1:size(u,2)
        u_TempFilt = transpose(squeeze(u(a, b, :)));
        u_freq = fft(u_TempFilt);
        
        % Condition to select equidistant points
        if mod(a, fix(size(u,1)/3)) == 0 && mod(b, fix(size(u,2)/3)) == 0 
            count = count + 1; % Increment the counter for subplots
            if count <= 9 % Limit to 9 subplots
                subplot(3, 3, count); % Create 3x3 subplots
                scale_x = [0 100 400 1200 2000];
                plot(freq_range, abs(u_freq(1:size(u,3)/2 + 1))./abs(max(u_freq)));
                xlabel('Frequency [Hz]');
                ylabel('Normalized Magnitude');
                title(['Point (', num2str(a), ',', num2str(b), ')']);
                set(gca, 'XTick', scale_x, 'XtickLabel', scale_x);
            end
        end
        
        % Extract the phase if needed
        if index ~= 0
            coef(a, b) = u_freq(index);
        else
            coef(a, b) = 0;
        end
    end
end

% Adjust figure size and layout if needed
set(fig, 'Position', [100, 100, 1200, 900]);


%% Creating 2D Bandpass filter and filtering
[f1,f2] = freqspace([size(u,1) size(u,2)],'meshgrid');
f1=f1*(2*pi*f_spacial_z); % Denormalizing the domain
f2=f2*(2*pi*f_spacial_x); % Denormalizing the domain
r = sqrt(f1.^2 + f2.^2);
Hd = ones([size(u,1) size(u,2)]); 
Hd((r< k_min*2)|(r>k_max*2)) = 0; % No need to normalize since the domain is in the same scale. Times '2' to get half the sampling frequency in the division
win=fspecial('gaussian',[size(u,1) size(u,2)],75); % Order 50 filter
win = win ./ max(win(:)); % Normalizing the window
h = fwind2(Hd,win); % Filter creation 

    % Ploting filter Window
win_size=101;
h_f=fftshift(fft2(h));
x_axis_sf=linspace(-fix(win_size/2),fix(win_size/2),win_size)*f_spacial_x*(10^-3)/size(u,2);
z_axis_sf=linspace(-fix(win_size/2),fix(win_size/2),win_size)*f_spacial_z*(10^-3)/size(u,1);
fig2=figure;
imagesc(x_axis_sf,z_axis_sf,abs(h_f(fix(size(u,1)/2)-fix(win_size/2):fix(size(u,1)/2)+fix(win_size/2),fix(size(u,2)/2)-fix(win_size/2):fix(size(u,2)/2)+fix(win_size/2)))); colormap jet; colorbar
xlabel('Spatial freq. [1/mm]');ylabel('Spatial freq. [1/mm]'); title('2-D FIR Bandpass filter');    
% saveas(fig2,'2DFIRFilter.png');
    
    %Ploting Input Signal Spectrum Window
plot_coef=fftshift(fft2(coef));
fig3=figure; 
imagesc(x_axis_sf,z_axis_sf,log10(abs(plot_coef(fix(size(u,1)/2)-fix(win_size/2):fix(size(u,1)/2)+fix(win_size/2),fix(size(u,2)/2)-fix(win_size/2):fix(size(u,2)/2)+fix(win_size/2))))); colormap jet; colorbar
xlabel('Spatial freq. [1/mm]');ylabel('Spatial freq. [1/mm]'); title('Spatial Frequency Spectrum of input signal');
% saveas(fig3,'Spatial_Frequency_Spectrum_input_signal.png');

    %Ploting Filter & Signal Profile Windows
x_value=fix(size(u,2)/2);y_value=fix(size(u,1)/2); window_l=50;
conv_factor_x=f_spacial_x*(10^-3);
conv_factor_z=f_spacial_z*(10^-3);
%Signal axial and lateral vectors
axial_value=fftshift(fft(coef(:,x_value)'));lateral_value=fftshift(fft(coef(y_value,:)));
axial_value=axial_value./max(axial_value(:));lateral_value=lateral_value./max(lateral_value(:)); % Normalizing
%Filter axial and lateral vectors
axial_valuef=fftshift(fft(h(:,x_value)'));lateral_valuef=fftshift(fft(h(y_value,:)));
axial_valuef=axial_valuef./max(axial_valuef(:));lateral_valuef=lateral_valuef./max(lateral_valuef(:)); % Normalizing
fig4=figure;
plot((-fix(window_l/2):fix(window_l/2))*conv_factor_z/size(u,1),abs(axial_value(y_value-fix(window_l/2):y_value+fix(window_l/2))),(-fix(window_l/2):fix(window_l/2))*conv_factor_z/size(u,1),abs(axial_valuef(y_value-fix(window_l/2):y_value+fix(window_l/2))));
xlabel('Spatial Freq.[1/mm]'); ylabel('Normalized Magnitude');title('Axial Profile');legend('data','filter')
% saveas(fig4,'AxialProfile.png');

fig5=figure;
plot((-fix(window_l/2):fix(window_l/2))*conv_factor_x/size(u,2),abs(lateral_value(x_value-fix(window_l/2):x_value+fix(window_l/2))),(-fix(window_l/2):fix(window_l/2))*conv_factor_x/size(u,2),abs(lateral_valuef(x_value-fix(window_l/2):x_value+fix(window_l/2))));
xlabel('Spatial Freq.[1/mm]'); ylabel('Normalized Magnitude');title('Lateral Profile');legend('data','filter')
% saveas(fig5,'LateralProfile.png');

    % Spacial Filtering
u_Filtered=filter2(h,coef);
%u_Filtered=imfilter(coef,h,'conv');
%u_EspFilt=fftshift(fft2(coef)).*Hd;
%u_Filtered=ifft2(u_EspFilt);

%% Ploting Original signal
    %Magnitude
fig6=figure;
imagesc(xx*10,yy*10,abs(coef)); colormap jet; colorbar
%imagesc(abs(coef)); colormap jet; colorbar
%conversion_x=(resol_x*10^3); % in mm
%conversion_y=(resol_z*10^3); % in mm
%addMMX=@(x) sprintf('%.1f',x*conversion_x);
%xticklabels(cellfun(addMMX,num2cell(xticks'),'UniformOutput',false));
%addMMY=@(x) sprintf('%.1f',x*conversion_y);
%yticklabels(cellfun(addMMY,num2cell(yticks'),'UniformOutput',false));
xlabel('Lateral [mm]');ylabel('Depth [mm]'); title('Original Signal Magnitude');
saveas(fig6,'OriginalSignalMagnitude.png');
    %Phase
fig7=figure;
imagesc(xx*10,yy*10,angle(coef)); colormap jet; colorbar
%imagesc(angle(coef)); colormap jet; colorbar
%conversion_x=(resol_x*10^3); % in mm
%conversion_y=(resol_z*10^3); % in mm
%addMMX=@(x) sprintf('%.1f',x*conversion_x);
%xticklabels(cellfun(addMMX,num2cell(xticks'),'UniformOutput',false));
%addMMY=@(x) sprintf('%.1f',x*conversion_y);
%yticklabels(cellfun(addMMY,num2cell(yticks'),'UniformOutput',false));
xlabel('Lateral [mm]');ylabel('Depth [mm]'); title('Original Signal Phase');
saveas(fig7,'OriginalSignalPhase.png');
%% Ploting filtered signal
    %Magnitude
fig8=figure;
imagesc(xx*10,yy*10,abs(u_Filtered)); colormap jet; colorbar
%conversion_x=(resol_x*10^3); % in mm
%conversion_y=(resol_z*10^3); % in mm
%addMMX=@(x) sprintf('%.1f',x*conversion_x);
%xticklabels(cellfun(addMMX,num2cell(xticks'),'UniformOutput',false));
%addMMY=@(x) sprintf('%.1f',x*conversion_y);
%yticklabels(cellfun(addMMY,num2cell(yticks'),'UniformOutput',false));
xlabel('Lateral [mm]');ylabel('Depth [mm]'); title('Filtered Signal Magnitude');
saveas(fig8,'FilteredSignalMagnitude.png');
    %Phase
fig9=figure;
imagesc(xx*10,yy*10,angle(u_Filtered)); colormap jet; colorbar
%conversion_x=(resol_x*10^3); % in mm
%conversion_y=(resol_z*10^3); % in mm
%addMMX=@(x) sprintf('%.1f',x*conversion_x);
%xticklabels(cellfun(addMMX,num2cell(xticks'),'UniformOutput',false));
%addMMY=@(x) sprintf('%.1f',x*conversion_y);
%yticklabels(cellfun(addMMY,num2cell(yticks'),'UniformOutput',false));
xlabel('Lateral [mm]');ylabel('Depth [mm]'); title('Filtered Signal Phase');
saveas(fig9,'FilteredSignalPhase.png');

    %Spatial Frequency domain
test=fftshift(fft2(u_Filtered));
x_axisplot=linspace(-fix(size(u,2)/2),fix(size(u,2)/2),size(u,2))*f_spacial_x*(10^-3)/size(u,2);
z_axisplot=linspace(-fix(size(u,1)/2),fix(size(u,1)/2),size(u,1))*f_spacial_z*(10^-3)/size(u,1);
fig10=figure; 
imagesc(x_axisplot,z_axisplot,log10(abs(test))); colormap jet;
xlabel('Spatial freq. [1/mm]');ylabel('Spatial freq. [1/mm]'); title('Filtered Signal Spectrum LogScale');    
saveas(fig10,'FilteredSignalSpectrum_LogScale.png');
%% Estimaging wave velocity

window_size=89; % must be an odd number
u1=size(u_Filtered);

z_pad=zeros([2*fix(window_size/2)+u1(1),2*fix(window_size/2)+u1(2)]);
z_pad((fix(window_size/2)+1):end-fix(window_size/2),(fix(window_size/2)+1):end-fix(window_size/2))=u_Filtered;

dz=3; delta_z=resol_z*dz;
dx=3; delta_x=resol_x*dx;
%d135=4; delta_135=resol*d135;
%d45=4; delta_45=resol*d45;
A_n=xcorr2(ones(window_size,window_size)); % For normalizing autocorrelation
wave_speed=zeros(u1(1),u1(2));

for a=1:u1(2)
    parfor b=1:u1(1)
        work_area=z_pad(b:window_size+b-1,a:window_size+a-1);
        auto_corr=xcorr2(work_area);
        auto_corr=auto_corr./A_n;
        
        z_vector=transpose(auto_corr(:,window_size));
        z_vector=z_vector./max(abs(z_vector)); % Normalizing
        k_z2=(10/((delta_z^2)*real(z_vector(window_size))))*(real(z_vector(window_size))-real(z_vector(window_size-dz)));
        
        x_vector=auto_corr(window_size,:);
        x_vector=x_vector./max(abs(x_vector)); % Normalizing
        k_x2=(15/(pi*(delta_x^2)*real(x_vector(window_size))))*(real(x_vector(window_size))-real(x_vector(window_size-dx)));
        
        %vector_45=transpose(diag(fliplr(auto_corr)));
        %vector_45=vector_45./max(abs(vector_45)); % Normalizing
        %k_45=(5/(sqrt(2*pi)*(delta_45^2)*real(vector_45(window_size))))*(real(vector_45(window_size))-real(vector_45(window_size-d45)));
        
        %vector_135=transpose(diag(auto_corr));
        %vector_135=vector_135./max(abs(vector_135)); % Normalizing
        %k_135=(5/(sqrt(2*pi)*(delta_135^2)*real(vector_135(window_size))))*(real(vector_135(window_size))-real(vector_135(window_size-d135)));
        
        %k_final=sqrt((k_z2+k_x2+k_45+k_135)/4);
        k_final=(sqrt(k_x2)+sqrt(k_z2))*0.5;
        if a==fix(size(u,2)/2) && b==fix(size(u,1)/2)
            if resol_x==resol_z
                fig11=figure;
                plot(resol_z*(10^3)*(-(window_size-1):window_size-1),real(z_vector),resol_x*(10^3)*(-(window_size-1):window_size-1),real(x_vector));
                grid on
                xlabel('size [mm]');ylabel('Normalized autocorrelation function');legend('axial','lateral'); title('Autocorrelation profiles');
                % saveas(fig11,'AutocorrelationProfiles.png');
            else
                fig11=figure;
                plot(resol_z*(10^3)*(-(window_size-1):window_size-1),real(z_vector));
                grid on
                xlabel('size [mm]');ylabel('Normalized autocorrelation function'); title('Autocorrelation Axial profile');
                % saveas(fig11,'AutocorrelationAxialProfile.png');
                
                fig12=figure;
                plot(resol_x*(10^3)*(-(window_size-1):window_size-1),real(x_vector));
                grid on
                xlabel('size [mm]');ylabel('Normalized autocorrelation function');title('Autocorrelation Lateral profile');
                % saveas(fig12,'AutocorrelationLateralProfile.png');
            end
            
            fig13=figure;
            x_scale_au=linspace(-window_size,window_size)*resol_x*(10^3); % in mm
            z_scale_au=linspace(-window_size,window_size)*resol_z*(10^3); % in mm
            imagesc(x_scale_au,z_scale_au,abs(auto_corr)); colormap jet; colorbar; 
            xlabel('Lateral [mm]');ylabel('Axial [mm]');title('Autocorrelation Map'); 
            % saveas(fig13,'AutocorrelationMap.png'); 
        end
        wave_speed(b,a)=(2*pi*f_v)/k_final;
    end
end

fig14=figure;
imagesc(xx*10,yy*10,abs(wave_speed)); colormap jet; colorbar;
caxis([0 5.5])
%conversion_x=(resol_x*10^3); % in pixel/mm
%conversion_y=(resol_z*10^3); % in pixel/mm
%addMMX=@(x) sprintf('%.1f',x*conversion_x);
%xticklabels(cellfun(addMMX,num2cell(xticks'),'UniformOutput',false));
%addMMY=@(x) sprintf('%.1f',x*conversion_y);
%yticklabels(cellfun(addMMY,num2cell(yticks'),'UniformOutput',false));
xlabel('Lateral [mm]');ylabel('Depth [mm]'); title('SWS Image');
% saveas(fig14,'SWS_image');
mean(mean(abs(wave_speed)))