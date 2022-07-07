% compute distances between neighbouring links in a time series chains
function dist = Chain_dist_neib(A,B, chain_idx, SubsequenceLength)
gt = length(A);
A = [A; B];
dist=[];
for i=1:length(chain_idx)-1
    cur = chain_idx(i);
    %disp(cur)
    next = chain_idx(i+1);
    %figure(figure2);
    %title([input_entity, ' ', var_name])
    %subplot(ceil(length(chain_idx)/3),3,i);
    curpattern=A(cur:(cur+SubsequenceLength-1));
    nextpattern = A(next:(next+SubsequenceLength-1));
    
    curpattern=zscore(A(cur:(cur+SubsequenceLength-1)),1);
    nextpattern=zscore(A(next:(next+SubsequenceLength-1)),1);
    dist(i)= norm(curpattern-nextpattern,2);
end

end