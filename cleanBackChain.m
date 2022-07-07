function chain_tmp = cleanBackChain(chain_tmp, A,B, SubseqLength, threshold)

chain_len = length(chain_tmp);
neib_dist = Chain_dist_neib(A, B, chain_tmp, SubseqLength);

mean_dist = cumsum(neib_dist)./(1:chain_len-1);

%[val, idx] = max(flip(neib_dist)>threshold);
[val, idx] = max(neib_dist>threshold);
if idx==1 & val==1
    chain_tmp=[];
elseif idx==1 & val==0
    return 
else
    chain_tmp = chain_tmp(1:idx);
end