function plotTSChain_normalized2(A,B, chain_idx, mp, mpi, gt, k,SubseqLength)
% the blue lines are chain A
% the purple lines are Chain B
% the blue dotted lines are the nearest neighbor in A for that node
gt = length(A);
A = [A; B];

% find chain in B k_th discord dist
chainB_mp = mp(chain_idx(chain_idx>gt));
[val_tmp,i_tmp] = sort(chainB_mp, 'desc');

dist_topk_discords = val_tmp(1:min(k,length(val_tmp)));
idx_topk = i_tmp(1:min(k,length(val_tmp)));
%disp(chainB_mp)
%disp(dist_topk_discords)
%disp(idx_topk)


for i=1:length(chain_idx)
    cur = chain_idx(i);
    %figure(figure2);
    %title([input_entity, ' ', var_name])
    subplot(4,ceil(length(chain_idx)/4), i);
    curpattern=zscore(A(cur:(cur+SubseqLength-1)),1);
    %curpattern=A(cur:(cur+SubseqLength-1));
    
    %plot(cur:(cur+SubseqLength-1),curpattern);
    
    % determine color of plot by incident
    if cur<gt
        plot(cur:(cur+SubseqLength-1),curpattern, 'b', 'LineWidth',2);
    else
        plot(cur:(cur+SubseqLength-1),curpattern, 'm', 'LineWidth',2);
           
            hold on
          
            curnn = mpi(cur);
            %disp(curnn)
            nnpattern =zscore(A(curnn:(curnn+SubseqLength-1)),1);
            plot(cur:(cur+SubseqLength-1),nnpattern, 'b:', 'LineWidth',2);
            hold off
        end
            
            
    end
    xlim([cur,cur+SubseqLength-1]);
    ylim([min(curpattern)-0.1 max(curpattern)+0.1]);
  
    %if i>1
        % hold on;
        % plot(cur:(cur+SubseqLength-1),-lastpattern+curpattern,'r');
   % end
   % figure(figure1);
   % plot(cur:(cur+SubseqLength-1),A(cur:(cur+SubseqLength-1)),'r');
end

