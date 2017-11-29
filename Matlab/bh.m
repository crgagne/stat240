function [ p_threshold ] = bh( ps, q )
%BH calculate p threshold for Benjamini-Hochberg FDR
%   ps: a list of p-values (doesn't have to be sorted)
%   q: desired FDR
%   keep all ps<=p_threshold. returns -1 if q cannot be met

ps = sort(ps,'ascend');
N = length(ps);
p_threshold = -1;
for i = 1:N
    if ps(i) <= i/N*q;
        p_threshold = ps(i);
    end
end 

end

