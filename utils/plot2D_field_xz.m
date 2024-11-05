figure
max_uz=max(uz_s(:));
min_uz=min(uz_s(:));
for ii=1:5:3000
%im=squeeze((uz_s(10,:,:,ii)-min_uz)/(max_uz-min_uz))';
im=squeeze((uz_s(10,:,:,ii)))';
imagesc(xx,zz,im);
xlabel('x (mm)'); ylabel('z (mm)');
axis image;
clim([-0.01 0.01]); %modify limits to see waves
title(['SWS - XZ plane - Frame ' num2str(ii)]);
grid on;
colormap('jet')
colorbar;
drawnow
pause(0.1)
end