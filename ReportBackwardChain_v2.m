function [Chain, ChainLength, ChainStart] = ReportBackwardChain_v2(A, B, MPindexLeft, SubseqLength, threshold)
% input: time series A, B, left MP index, noise thereshold
% output: chain cell, chain length, position of chain start


MatrixProfileLength = length(MPindexLeft);

ChainLength = zeros(MatrixProfileLength,1);
Chain  = cell(length(MPindexLeft),1);


for i = MatrixProfileLength:-1:SubseqLength

        cur=i;
        count=1;
        Chain_tmp = i;
        while MPindexLeft(cur)>0 && (~isnan(MPindexLeft(cur)))
            cur=MPindexLeft(cur);
            Chain_tmp = [Chain_tmp cur];
            count=count+1;
        end
        ChainStart(i) = cur;
        %ChainLength(i)=count;
        Chain{i} = Chain_tmp;
end


for i = MatrixProfileLength:-1:SubseqLength
        Chain{i}=fliplr(cleanBackChain(Chain{i}, A,B, SubseqLength, threshold));
        ChainLength(i)=length(Chain{i});
        if isempty(Chain{i})
            ChainStart(i) = 0;
        else
            ChainStart(i) = Chain{i}(end);
        end
end

end