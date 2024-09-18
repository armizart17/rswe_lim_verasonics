function [cx,cz,c,VelF] = elastoTRE(Vel,dinf,lambda_cut,type)
% function [cx,cz,c,VelF] = elastoTRE(Vel,dinf,lambda_cut,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the shear wave speed of a region with the 
% UDELAR method. 
% 
% Inputs:  
%          Vel         - particle velocity
%          dinf        - for resolution, should include in [m]
%                   dinf.dx (Lateral)
%                   dinf.dz (Axial)
%          lambda_cut  - wavelength cut for spatial filtering
%          type        - kind of filter for filt2, i.e. 'lp', 'hp', 'bp', 'bs'
%
% Outputs:  
%          cx          - speed lateral
%          cx          - speed axial
%          c           - speed
%          VelF        - particle velocity after filtering
% Author: E.A.M.Z based on LIM Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
dx = dinf.dx;
dz = dinf.dz;
res = 1000*(dx+dz)/2; % promedio entre resoluciones segun z y segun x
% lambda_cut = 3; % parámetro de filtro para la función filt2
% type = 'lp';
nptos = 1;
% nptos = 5; emz
disp("a")
if nargin<5
    VelF = zeros(size(Vel));
    for ii=1:size(Vel,3)
         VelF(:,:,ii) = filt2(Vel(:,:,ii),res,lambda_cut,type);
    %       VelF(:,:,ii)=imgaussfilt(squeeze(Vel(:,:,ii)),1.5); 
%         VelF(:,:,ii)=medfilt2(Vel(:,:,ii),[15 15]);
    end
end
% for frames = 1:250
%     subplot(121),imagesc(Vel(:,:,frames))
%     subplot(122),imagesc(VelF(:,:,frames))
%     pause;
% end

% [nz,nx,nt] = size(VelF);

%% TRE ancho focal

% Dpl=VelF(1:4:end,:,:);
% EDpl=1./sqrt(squeeze(sum(Dpl.*Dpl,3)));
% % [nx nz nt]=size(Dpl);
% % To=ceil(nt/2)+1;   % calculo del maximo de la focalizacion temporal
% 
% Elasto90=zeros(size(Dpl,1),size(Dpl,2));
% Elasto80=zeros(size(Dpl,1),size(Dpl,2));
% %RTij=zeros(size(Dpl,1),size(Dpl,2));
% for io=1:size(Dpl,1);
%     for jo=1:size(Dpl,2)
% %         io=15; % zeta
% %         jo=32; % equis
% 
%         Dpl_ij=repmat(Dpl(io,jo,:),size(Dpl,1),size(Dpl,2));
%         RTij=sum(Dpl.*Dpl_ij,3).*EDpl*EDpl(io,jo);
%         Elasto90(io,jo)=length(find(RTij>0.90));
%         Elasto80(io,jo)=length(find(RTij>0.80));
% %         jo
%     end
%     io
% end
% 


%% TRE a la Stefan (gradient)

% %1) FILTRADO ESPACIAL DE DATOS
Dpl = VelF;
%Dpl=cumsum(Velocity,3);
%Dpl=u;


% VelF=Vel;
[ez,ex,vpart] = gradient(Dpl,nptos*dx,nptos*dz,1/dinf.PRFe); % caluclo de gradiente
% kernel=1;
% volk=3*3;
% ez=convn(ez,ones(kernel,kernel,1)/volk);
% ex=convn(ex,ones(kernel,kernel,1)/volk);
% vpart=convn(vpart,ones(kernel,kernel,1)/volk);
% Dpl=convn(Dpl,ones(kernel,kernel,1)/volk);

% for ii=1:size(Vel,3)
%     ez(:,:,ii)=imgaussfilt(squeeze(ez(:,:,ii)),2);
%     ex(:,:,ii)=imgaussfilt(squeeze(ex(:,:,ii)),2);
% %     ez(:,:,ii)=filt2(ez(:,:,ii),res,lambda_cut,type);
% %     ex(:,:,ii)=filt2(ex(:,:,ii),res,lambda_cut,type);
% end


% for ii=11:size(ez,3);
% %     ez2(:,:,ii)=conv2(squeeze(ez(:,:,ii)),ones(4,4)/16);
%     ez(:,:,ii)=medfilt2(ez(:,:,ii),[5 5]);
%     ex(:,:,ii)=medfilt2(ex(:,:,ii),[5 5]);
% end


% RTij = sum(Dpl.*Dpl,3);
ez_TR = sum(ez.*ez,3);
ex_TR = sum(ex.*ex,3);
vpart_TR = sum(vpart.*vpart,3);

% clear ex ez vpart

cz = sqrt(vpart_TR./ez_TR);
cx = sqrt(vpart_TR./ex_TR);
c = sqrt(vpart_TR./(ez_TR+ex_TR));
% lambda_x = 2*pi*sqrt(RTij./ex_TR);
% lambda_z = 2*pi*sqrt(RTij./ez_TR);
% lambda = 2*pi*sqrt(RTij./(ex_TR+ez_TR));


% figure; imagesc((medfilt2(cx,[5 5]))); colormap jet