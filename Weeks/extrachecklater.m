load ('P:\emiranda\proj\rswe\kwave_v1\simRSWF_2D\RSWF_2D_v1.mat')

ux = gather(sensor_data.ux);
uz = gather(sensor_data.uy);

figure, 
subplot(121)
    imagesc(abs(ux))
    title('ux')
subplot(122)
    imagesc(abs(uz))
    title('uz')

    %%

uz_s = sensor_data.uz_split_s;

% modify according to the size of the sensor (Nx,..,Nz,frames);
uz_s = reshape(uz_s,2*thick-1,Ny,Nz,[]);


%%

uy_s = gather(sensor_data.uy_split_s);
Nx = 80;
Ny = 80;
% modify according to the size of the sensor (Nx,..,Nz,frames);
uz_s = reshape(uy_s, Nx, Ny,[]);

%%
load ('P:\emiranda\proj\rswe\kwave_v1\simRSWF\RSWF_v1p6_homo.mat')