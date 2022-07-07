function [Chain, Chain_score, L] = TSC_join(A,B, SubseqLength, f, run_mp)

A(isnan(A))=0;
B(isnan(B))=0;
gt = length(A);

f1 = [f '_mp'];


if run_mp
    disp('compute matrix profiles.....')
    % compute AB join matrix profile
    [mp, mpi] = MatrixProfileAB([A;B], SubseqLength,gt);
    
    mpia = mpi(1:length(A))-length(A);
    mpib = mpi(length(A)+1:end);
    mpia(mpia<1)=0;
    mpia(length(A)-SubseqLength+1:length(A))=-1;
    % compute left and right matrix profiles for A and B
    [MPLeft_A, MPRight_A, MPindexLeft_A, MPindexRight_A]= ComputeMP(A, SubseqLength);
    [MPLeft_B, MPRight_B, MPindexLeft_B, MPindexRight_B]= ComputeMP(B, SubseqLength);
    disp(['finish computing matrix profiles and saved in' f1 '.mat'])
    save(f)
else
    load(f1)
    disp(['loaded matrix profiles of ' f1 '.mat'])
end

% case study parameters
qb = 0.5; % noise quantile in TB
qa = 0.5; % noise quantile in TA
%theta_AB = 20; % AB joint node max distance
theta_A = quantile(MPLeft_A, qa);
theta_B = quantile(MPRight_B, qb);
k=3; % top k Chain B in ranking

%theta_AB = 1000;

% compute forward chain in TB
[ChainB, ChainLenB, ChainEndB] = ReportForwardChain_v2([], B, MPindexRight_B, SubseqLength, theta_B);
% backward chain in TA
[ChainA, ChainLenA, ChainStartA] = ReportBackwardChain_v2(A, [], MPindexLeft_A, SubseqLength, theta_A);
% combine chain
[Chain,ChainAc, L] = CombineABChain(ChainA,ChainB,mpia,mpib,mp,length(A));

% ranking score
Chain_score = ComputeChainRank(Chain, ChainAc, mp, A, B, length(A), SubseqLength, k);

% max chain
[a,b] = max(Chain_score);

% compute precursor
topChain = Chain{b};
topChain = topChain(topChain>length(A));

top_neib_dist = Chain_dist_neib(A,B, topChain, SubseqLength);
chainB_mp = mp(topChain);
disp(['Chain detected:' num2str(Chain{b})]) 
disp('Done.');

figure

plotTSChain_normalized2(A,B, Chain{b}, mp, mpi, length(A), k,SubseqLength)
