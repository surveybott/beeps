function [smoothedSignal] = sliding_window(in,window)

for i = 1:length(in)
    
    if length(in(i:end)) > window
    
        smoothedSignal(i,1) = nanmean(in(i:i+window));
    
    else
        
        smoothedSignal(i,1) = nanmean(in(i:end));
        
    end


end

