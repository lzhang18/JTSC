function Demo_JTSC()
load('demo_data')
%the data is in demo_data.mat
%A: real-value array of normal behavior data
%B: real-value array of poptentially abnormal data 
%f: file name for storing matrix profile (if data is 'DataAB', then f='DataAB')
%run_mp: 1 if we need to rerun matrix profile; 
%        0 if we stored matrix profile before
%SubseqLength: subsequence length of time series chain. Default value is set as 45.  

% noise cut parameters for TSC_join is in Line 45-49 TSC_Join

[Chain, Chain_score, L] = TSC_join(A,B, SubseqLength, f, 1);