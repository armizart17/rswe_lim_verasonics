%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functition that finds the position in which the closest value to the
% vibration frequency is located in the frequency range.
% 
% Inputs:  f_v         - vibration freqeuncy
%          freq_range  - The frequency range vector 
%
% Outputs: index       - position in the frequency range vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function index = find_position(f_v,freq_range)
    index_0 = min(abs(f_v-freq_range)); 
    index = find(freq_range == f_v + index_0);
    if isempty(index)
        index=find(freq_range==f_v-index_0);
    end
end