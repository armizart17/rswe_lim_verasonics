function [A, P, mat_complex] = amp_phs(u, ct, ct2)
% function [A, P, mat_comp] = amp_phs(u, ct, ct2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnitude and phase extraction from temporal signals
% Author: Edited by EMZ from LIM repository
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A = zeros(size(u,1),size(u,2));
    P = zeros(size(u,1),size(u,2));
    for ii = 1:size(u,1)
        for jj = 1:size(u,2)
            s1 = squeeze(u(ii,jj,:));
            F1 = fftshift(fft(s1,ct2)); 
    %         M = abs(F1);             % magnitude
            P(ii,jj) = angle(F1(ct));
            A(ii,jj) = abs(F1(ct));
        end
    end
    mat_complex = A.*exp(1i*P);
end 