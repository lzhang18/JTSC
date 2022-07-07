function [diviance_score, dist_topk_discord] = ComputeChainScore(chain_idx, mp, A, B, gt, SubsequenceLength, k)

chain_idx_A = chain_idx(chain_idx < gt);
chain_idx_B = chain_idx(chain_idx >= gt);

% find chain in B k_th discord dist
chainB_mp = mp(chain_idx_B);
[val_tmp,i_tmp] = sort(chainB_mp, 'desc');
dist_topk_discord = val_tmp(min(k,length(val_tmp)));


dist_neib_A = Chain_dist_neib(A,B, chain_idx_A, SubsequenceLength);

dist_node_join = Chain_dist_neib(A, B, [chain_idx_A(end) chain_idx_B(1)], SubsequenceLength);


dist_dev_B = Chain_dist_neib(A,B, [chain_idx_B(1) chain_idx_B(end)], SubsequenceLength);


kth_discord_score = (dist_topk_discord+1)/(max([dist_neib_A, dist_node_join])+1);

diviance_score = kth_discord_score* dist_dev_B; 

end