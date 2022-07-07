function [Chain, ChainLength, ChainEnd] = ReportForwardChain_v2(A, B, MPindexRight_B, SubseqLength, threshold)
% input: time series A, B, right MP index, noise thereshold
% output: chain cell, chain length, position of chain end

MatrixProfileLength = length(MPindexRight_B);

ChainLength = zeros(MatrixProfileLength,1);
Chain  = cell(length(MPindexRight_B),1);


for i = 1:MatrixProfileLength-SubseqLength+1
        cur=i;
        count=1;
        Chain_tmp = i;
        
        while MPindexRight_B(cur)>0 && (~isnan(MPindexRight_B(cur)))
            cur=MPindexRight_B(cur);
            Chain_tmp = [Chain_tmp cur];
            count=count+1;
        end
        ChainEnd(i) = cur;
        Chain{i} = Chain_tmp;
        ChainLength(i)=count;
            

end

for i = 1:MatrixProfileLength-SubseqLength+1
        Chain{i}=cleanForwardChain(Chain{i}, A,B, SubseqLength, threshold);
        ChainLength(i)=length(Chain{i});
        if isempty(Chain{i})
            ChainEnd(i) = 0;
        else    
            ChainEnd(i) = Chain{i}(end);
        end
end