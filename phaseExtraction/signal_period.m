function u_new = signal_period(fv, PRFe, u)
% function u_new = signal_period(fv, PRFe, u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After particle velocity estimation, a integer amount of periods are used,in
% in this case we take 10 periods of vibration signal
% Inputs:  
%           fv      : vibration freqeuncy
%           PRFe    : pulse repetition frequency         
%           u       : particle velocity data
% Outputs: 
%           u_new   : particle velocity data with 10 periods
% Author. Edited by EMZ from LIM repository
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P = 1/fv;
    pm = 1/PRFe;
    limit = round(size(u,3)/(PRFe/fv));
    N = limit-1; %10 
    fpoint = round(P*N/pm);
    spoint = 10;    %40
    
    while (spoint+fpoint-1) > size(u,3)
        N = N-1;
        fpoint = round(P*N/pm);     
    end
    
        u_new = u(:,:,spoint:spoint+fpoint-1);
    % u_new=u(:,:,spoint:end-1);
end