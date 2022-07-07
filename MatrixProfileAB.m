% The prototype for interactive matrix profile calculation
% Chin-Chia Michael Yeh / Eamonn Keogh 01/26/2016
%
% [matrixProfile, profileIndex, motifIndex, discordIndex] = ...
%     interactiveMatrixProfile(data, subLen, split);
% Output:
%     matrixProfile: matrix profile of the input data (vector)
%     profileIndex: matrix profile index of the input data (vecto r)
%     motifIndex: index of the first, second, and third motifs and their
%                 associated neighbors when stopped (3x2 cell)
%                 +-----------------------+------------------------+
%                 | indices for 1st motif | neighbors of 1st motif |
%                 +-----------------------+------------------------+
%                 | indices for 2nd motif | neighbors of 2nd motif |
%                 +-----------------------+------------------------+
%                 | indices for 3rd motif | neighbors of 3rd motif |
%                 +-----------------------+------------------------+
%     discordIndex: index of discords when stopped (vector)
% Input:
%     data: input time series (vector)
%     subLen: subsequence length (scalar)
%     split: the joins must cross this index, put inf to remove this
%            constrain (scalar)
%
% Chin-Chia Michael Yeh, Yan Zhu, Liudmila Ulanova, Nurjahan Begum,
% Yifei Ding, Hoang Anh Dau, Diego Furtado Silva, Abdullah Mueen, and
% Eamonn Keogh, "Matrix Profile I: All Pairs Similarity Joins for Time
% Series: A Unifying View that Includes Motifs, Discords and Shapelets,"
% ICDM 2016, http://www.cs.ucr.edu/~eamonn/MatrixProfile.html
%

function [matrixProfile, profileIndex] = ...
    MatrixProfileAB(data, subLen, split)
%% options for the algorithm
excZoneLen = round(subLen * 0.5);

anytimeMode = 2; % 1: original with mass O(n^2 log n);
% 2: new diagnal method O(n^2)
dataLen = length(data);

%% locate nan and inf
proLen = dataLen - subLen + 1;
isSkip = false(proLen, 1);
for i = 1:proLen
    if any(isnan(data(i:i + subLen - 1))) || ...
            any(isinf(data(i:i + subLen - 1)))
        isSkip(i) = true;
    end
end
data(isnan(data) | isinf(data)) = 0;

%% preprocess for matrix profile
[dataFreq, dataMu, dataSig] = massPre(data, dataLen, subLen);
matrixProfile = inf(proLen, 1);
profileIndex = zeros(proLen, 1);
if anytimeMode == 1
    idxOrder = randperm(proLen);
elseif anytimeMode == 2
    idxOrder = excZoneLen + 1:proLen;
    idxOrder = idxOrder(randperm(length(idxOrder)));
end

%% main loop

for i = 1:length(idxOrder)
    idx = idxOrder(i);
    if isSkip(idx)
        continue
    end
    drawnow;
    
    % compute the distance profile and update matrix profile
    query = data(idx:idx+subLen-1);
    if anytimeMode == 1
        distProfile = mass(dataFreq, query, dataLen, subLen, ...
            dataMu, dataSig, dataMu(idx), dataSig(idx));
        distProfile = abs(distProfile);
        distProfile = sqrt(distProfile);
        
        distProfile(isSkip) = inf;
        excZoneStart = max(1, idx - excZoneLen);
        excZoneEnd = min(proLen, idx + excZoneLen);
        distProfile(excZoneStart:excZoneEnd) = inf;
        
        updatePos = distProfile < matrixProfile;
        profileIndex(updatePos) = idx;
        matrixProfile(updatePos) = distProfile(updatePos);
        [matrixProfile(idx), profileIndex(idx)] = min(distProfile);
    elseif anytimeMode == 2
        distProfile = diagonalDist(...
            data, idx, dataLen, subLen, proLen, dataMu, dataSig);
        distProfile = abs(distProfile);
        distProfile = sqrt(distProfile);
        
        pos1 = idx:proLen;
        pos2 = 1:proLen - idx + 1;
        
        if ~isinf(split)
            distProfile = distProfile(pos2 <= split & pos1 > split);
            pos1Split = pos1(pos2 <= split & pos1 > split);
            pos2Split = pos2(pos2 <= split & pos1 > split);
            pos1 = pos1Split;
            pos2 = pos2Split;
        end
        
        updatePos = matrixProfile(pos1) > distProfile;
        profileIndex(pos1(updatePos)) = pos2(updatePos);
        matrixProfile(pos1(updatePos)) = distProfile(updatePos);
        updatePos = matrixProfile(pos2) > distProfile;
        profileIndex(pos2(updatePos)) = pos1(updatePos);
        matrixProfile(pos2(updatePos)) = distProfile(updatePos);
        
        matrixProfile(isSkip) = inf;
        profileIndex(isSkip) = 0;
    end
    
    if(mod(i,1000)==0)
        disp(i);
    end
    
end


function [motifIdxs, matrixProfileCur] = findMotifs(...
    matrixProfileCur, profileIndex, dataLen, subLen, proLen, ...
    data, dataFreq, dataMu, dataSig, isSkip, excZoneLen, radius)
motifIdxs = cell(3, 2);
for i = 1:3
    [motifDistance, minIdx] = min(matrixProfileCur);
    motifDistance = motifDistance ^ 2;
    motifIdxs{i, 1} = sort([minIdx, profileIndex(minIdx)]);
    motifIdx = motifIdxs{i, 1}(1);
    query = data(motifIdx:motifIdx + subLen - 1);
    
    distProfile = mass(dataFreq, query, ...
        dataLen, subLen, dataMu, dataSig, ...
        dataMu(motifIdx), dataSig(motifIdx));
    distProfile = abs(distProfile);
    distProfile(distProfile > motifDistance * radius) = inf;
    motifZoneStart = max(1, motifIdx - excZoneLen);
    motifZoneEnd = min(proLen, motifIdx + excZoneLen);
    distProfile(motifZoneStart:motifZoneEnd) = inf;
    motifIdx = motifIdxs{i, 1}(2);
    motifZoneStart = max(1, motifIdx - excZoneLen);
    motifZoneEnd = min(proLen, motifIdx + excZoneLen);
    distProfile(motifZoneStart:motifZoneEnd) = inf;
    distProfile(isSkip) = inf;
    [distanceOrder, distanceIdxOrder] = sort(distProfile, 'ascend');
    motifNeighbor = zeros(1, 10);
    for j = 1:10
        if isinf(distanceOrder(1)) || length(distanceOrder) < j
            break;
        end
        motifNeighbor(j) = distanceIdxOrder(1);
        distanceOrder(1) = [];
        distanceIdxOrder(1) = [];
        distanceOrder(abs(distanceIdxOrder - motifNeighbor(j)) < ...
            excZoneLen) = [];
        distanceIdxOrder(abs(distanceIdxOrder - motifNeighbor(j)) < ...
            excZoneLen) = [];
    end
    motifNeighbor(motifNeighbor == 0) = [];
    motifIdxs{i, 2} = motifNeighbor;
    
    removeIdx = cell2mat(motifIdxs(i, :));
    for j = 1:length(removeIdx)
        removeZoneStart = max(1, removeIdx(j) - excZoneLen);
        removeZoneEnd = min(proLen, removeIdx(j) + excZoneLen);
        matrixProfileCur(removeZoneStart:removeZoneEnd) = inf;
    end
end


function pushDiscardBtn(src, ~, btnNum)
mainWindowFig = get(src, 'parent');
mainWindow = get(mainWindowFig, 'userdata');
mainWindow.discardIdx = [mainWindow.discardIdx, ...
    mainWindow.motifIdxs{btnNum, 1}];
for i = 1:3
    set(mainWindow.discardBtn(i), 'enable', 'off');
end
set(mainWindow.fig, 'userdata', mainWindow);


function pushStopBtn(src, ~)
mainWindowFig = get(src, 'parent');
mainWindow = get(mainWindowFig, 'userdata');
mainWindow.stopping = true;
for i = 1:3
    set(mainWindow.discardBtn(i), 'enable', 'off');
end
set(src, 'enable', 'off');
set(mainWindow.fig, 'userdata', mainWindow);


% The following two functions are modified from the code provided in the
% following URL
% http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [dataFreq, dataMu, dataSig] = massPre(data, dataLen, subLen)
data(dataLen + 1:(subLen + dataLen)) = 0;
dataFreq = fft(data);
dataCumsum = cumsum(data);
data2Cumsum =  cumsum(data .^ 2);
data2Sum = data2Cumsum(subLen:dataLen) - ...
    [0; data2Cumsum(1:dataLen - subLen)];
dataSum = dataCumsum(subLen:dataLen) - ...
    [0; dataCumsum(1:dataLen - subLen)];
dataMu = dataSum ./ subLen;
data2Sig = (data2Sum ./ subLen) - (dataMu .^ 2);
dataSig = sqrt(data2Sig);


function distProfile = mass(dataFreq, query, ...
    dataLen, subLen, dataMu, dataSig, queryMu, querySig)
query = query(end:-1:1);
query(subLen+1:(subLen+dataLen)) = 0;
queryFreq = fft(query);
productFreq = dataFreq .* queryFreq;
product = ifft(productFreq);
distProfile = 2 * (subLen - ...
    (product(subLen:dataLen) - subLen * dataMu * queryMu) ./ ...
    (dataSig * querySig));


function distProfile = diagonalDist(...
    data, idx, dataLen, subLen, proLen, dataMu, dataSig)
xTerm = ones(proLen - idx + 1, 1) * ...
    (data(idx:idx + subLen - 1)' * data(1:subLen));
mTerm = data(idx:proLen - 1) .* ...
    data(1:proLen - idx);
aTerm = data(idx + subLen:end) .* ...
    data(subLen + 1:dataLen - idx + 1);
if proLen ~= idx
    xTerm(2:end) = xTerm(2:end) - cumsum(mTerm) + cumsum(aTerm);
end

distProfile = (xTerm - ...
    subLen .* dataMu(idx:end) .* dataMu(1:proLen - idx + 1)) ./ ...
    (subLen .* dataSig(idx:end) .* dataSig(1:proLen - idx + 1));
distProfile = 2 * subLen * (1 - distProfile);


function x = zeroOneNorm(x)
x = x - min(x(~isinf(x) & ~isnan(x)));
x = x / max(x(~isinf(x) & ~isnan(x)));


function mainResize(src, ~)
mainWindow = get(src, 'userdata');
figPosition = get(mainWindow.fig, 'position');
axGap = 38;
axesHeight = round((figPosition(4) - axGap * 5 - 60) / 6);
set(mainWindow.dataAx, 'position', ...
    [30, 5 * axesHeight+5 * axGap + 30, figPosition(3) - 60, axesHeight]);
set(mainWindow.profileAx, 'position', ...
    [30, 4 * axesHeight+4 * axGap + 30, figPosition(3) - 60, axesHeight]);
set(mainWindow.discordAx, 'position', ...
    [30, 30, figPosition(3) - 160, axesHeight]);
set(mainWindow.stopBtn, 'position', ...
    [figPosition(3) - 120, 30, 90, 20]);
set(mainWindow.dataText, 'position', ...
    [30, 6 * axesHeight + 5 * axGap + 30, figPosition(3) - 60, 18]);
set(mainWindow.profileText, 'position', ...
    [30, 5 * axesHeight + 4 * axGap + 30, figPosition(3) - 60, 18]);
set(mainWindow.discordText, 'position', ...
    [30, 1 * axesHeight + 30, figPosition(3) - 160, 18]);
for i = 1:3
    set(mainWindow.motifAx(i), 'position', ...
        [30, (4 - i) * axesHeight + (4 - i) * axGap + 30, ...
        figPosition(3) - 160, axesHeight]);
end
for i = 1:3
    set(mainWindow.motifText(i), 'position', ...
        [30, (5 - i) * axesHeight + (4 - i) * axGap + 30, ...
        figPosition(3) - 160, 18]);
end
for i = 1:3
    set(mainWindow.discardBtn(i), 'position', ...
        [figPosition(3) - 120, ...
        (4 - i) * axesHeight + (4 - i) * axGap + 30, 90, 20]);
end