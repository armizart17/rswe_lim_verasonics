res_y = dz;
res_x = dx;

% sub_vz = vz(ev_al_pos_z(200) + search_area_z, ev_al_pos_x(90) + search_area_x); % JO DATA
sub_vz = vz(ev_al_pos_z(200) + search_area_z, ev_al_pos_x(90) + search_area_x);% LIM DATA

correc = xcorr2(ones(win(1), win(2)));

M_filtered = sub_vz;
corre = xcorr2(M_filtered);

corre = corre./correc;
    
corre = corre./(abs(corre(round((size(corre,1)+1)/2),round((size(corre,2)+1)/2))));
y_v = (-floor(size(corre,1)/2):1:floor(size(corre,1)/2))*res_y;
x_v = (-floor(size(corre,2)/2):1:floor(size(corre,2)/2))*res_x;

aux_rot = zeros(size(corre));
for angle = 1:1:360
    aux_rot = aux_rot+imrotate(real(corre),angle,"nearest","crop");
    sim_autocorr = aux_rot/360;
end

R_cc_ax = real(corre(:,round((size(corre,2)+1)/2)));
I_cc_ax = imag(corre(:,round((size(corre,2)+1)/2)));
R_cc_la = real(corre(round((size(corre,1)+1)/2),:));
I_cc_la = imag(corre(round((size(corre,1)+1)/2),:));

R_cc_1d = sim_autocorr(round((size(corre,1)+1)/2),:);

figure, imagesc(real(sub_vz))
title('Window')

figure, 
subplot(1,3,1)
plot(y_v, R_cc_ax), title('Axial'), grid minor;

subplot(1,3,2)
plot(x_v, R_cc_la), title('Lateral'), grid minor;


subplot(1,3,3)
plot(x_v,R_cc_1d), title('Sim360'), grid minor;



%%
track = 0.004;
    location1 = (abs(y_v+track)==min(abs(y_v+track)));  
    loc1(1) = find(location1==1);
    location1 = (abs(y_v-track)==min(abs(y_v-track)));  
    loc1(2) = find(location1==1);
    
    y_v2 = y_v(loc1(1):loc1(2));
    R_cc_1 = R_cc_ax(loc1(1):loc1(2));
    I_cc_1 = I_cc_ax(loc1(1):loc1(2));
    % [~, idz1] = min(abs(real(R_cc_1((length(R_cc_1)+1)*0.5:end))-0.9));
    % [~, idz2] = min(abs(real(R_cc_1(0<y_v2 & y_v2<4e-3))-0.15));
    % limit1z=y_v2((length(R_cc_1)+1)*0.5 + idz1-1);
    % limit2z=y_v2((length(R_cc_1)+1)*0.5 + idz2-1);

track = 0.003;    
    location1 = (abs(x_v+track)==min(abs(x_v+track)));  
    loc1(1) = find(location1==1);
    location1 = (abs(x_v-track)==min(abs(x_v-track)));  
    loc1(2) = find(location1==1);
    
    x_v2 = x_v(loc1(1):loc1(2));
    R_cc_2 = R_cc_la(loc1(1):loc1(2));
    I_cc_2 = I_cc_la(loc1(1):loc1(2));

figure, 
subplot(1,3,1)
plot(y_v2,R_cc_1), title('Cut Axial'), grid minor

subplot(1,3,2)
plot(x_v2,R_cc_2), title('Cut Lateral'), grid minor

subplot(1,3,3)
plot(x_v,R_cc_1d), title('Sim360'), grid minor


    
    % [~, idx1] = min(abs(real(R_cc_2((length(R_cc_2)+1)*0.5:end))-0.9));
    % [~, idx2] = min(abs(real(R_cc_2(0<x_v2 & x_v2<4e-3))-0.15));
    % limit1x=x_v2((length(R_cc_2)+1)*0.5 + idx1-1);
    % limit2x=x_v2((length(R_cc_2)+1)*0.5 + idx2-1);
    
    [fitresult1, gof1] = createFit_atenuation_axial_RSWE(y_v2, R_cc_1',I_cc_1',-0.01,0.01);
    [fitresult2, gof2] = createFit_atenuation_lateral_RSWE(x_v2, R_cc_2,I_cc_2,-0.006,0.006);
    [fitresult3, gof3] = createFit_atenuation_1d_RSWE (x_v, R_cc_1d);
    
    R_squeare_real_axial = gof1(1).rsquare;
    R_squeare_real_later = gof2(1).rsquare;
    
    k_axial = fitresult1.k;
    k_lateral = fitresult2.k;
    
    R_squeare_real_1d = gof3(1).rsquare;
    k_1d = fitresult3.k;