function [Chain_score, discord_score]= ComputeChainRank(Chain, ChainAc, mp, A, B, gt, SubsequenceLength, k)
Chain_score = zeros([length(Chain), 1]);
discord_score = zeros([length(Chain), 1]);
for i=1:length(Chain)
    disp(i)
    chain_min = min(Chain{i});
    chain_max = max(Chain{i});
    if isempty(ChainAc{i}) 
        Chain_score(i)=-1;
        continue
    elseif chain_min>gt || chain_max<gt
        Chain_score(i)=-1;
        continue
    else
        [Chain_score(i), discord_score(i)] = ComputeChainScore(Chain{i}, mp, A, B, gt, SubsequenceLength, k);
    end
end

end