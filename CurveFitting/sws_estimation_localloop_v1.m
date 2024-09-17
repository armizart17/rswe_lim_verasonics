function [Raxial,Rlateral,k_axial,k_lateral] = sws_estimation_localloop_v1(s_2D,dx,dy,N,n,correc)
% function [Raxial,Rlateral,k_axial,k_lateral] = sws_estimation_localloop_v1(s_2D,dx,dy,N,n,correc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D cross corralation procedure and second moment equation application

    al_pos_x = round(N/2):n-round(N/2); % allow locations to use in lateral direction
    ev_index_x = 1:floor(length(al_pos_x)); % even index vector 
    ev_al_pos_x = al_pos_x(ev_index_x); % even positions from allow locations in lateral direction
    search_area_x = -round(N/2)+1:round(N/2)-1; % search index from kernel:lateral

    Raxial = zeros(1,length(ev_al_pos_x));
    Rlateral = zeros(1,length(ev_al_pos_x));
    k_axial = zeros(1,length(ev_al_pos_x));
    k_lateral = zeros(1,length(ev_al_pos_x));
    mode = 1; %1 for xcorr, 2 for w-k theorem

    % for l = 1:length(ev_al_pos_x)
    parfor l = 1:length(ev_al_pos_x)
        [Raxial(:,l),Rlateral(:,l), k_axial(:,l),k_lateral(:,l)] = test_reverb_cf(s_2D(:,ev_al_pos_x(l)+search_area_x),dx,dy,mode,correc);
    end
end

