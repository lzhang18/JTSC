function [Chain, ChainAc, L] = CombineABChain(ChainA,ChainB,MPia,MPib,mp,offset, threshold)
% Combine chain
% input: 
% ChainA,ChainB --chainA, B, 
% MPia,MPib,mp -- matrix profiles
% offset -- mp overlapping
% threshold -- max allowed distance between first and last chains.
%              usually not useful unless for extreme outliers. Default=20.
% output: 
% Chain -- all combined chain in cell
% ChainAc -- Chain in A
% L -- length of Chain

if nargin==6
    threshold = 20;
end

Chain = cell(length(ChainB),1);
ChainAc = cell(length(ChainB),1);
L=zeros(length(ChainB),1);
for i = 1:length(ChainB)
    
    if(MPia(MPib(i))==i && mp(i+offset)<threshold)
        end_pos = MPib(i);
        disp(i)
        Chain{i}=[ChainA{end_pos} ChainB{i}+offset];
        ChainAc{i}=ChainA{end_pos};
        L(i) = length(Chain{i});
    else
        continue
    end
end