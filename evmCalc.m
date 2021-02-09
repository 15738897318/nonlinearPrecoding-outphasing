function evmOut = evmCalc(sF,hT,aOut,numUsedsubcarrier,usedSubcarrierIdx)
% - Calculate EVM given channel, payload data, and precoder output - %

numTrial = numel(sF);
N        = size(aOut{1},2);

if ~exist('usedSubcarrierIdx','var')
    usedSubcarrierIdx   = (-numUsedsubcarrier/2:numUsedsubcarrier/2-1) + N/2;
end
    
for ii = 1:numTrial
    hF = fft(hT{ii},N,3);
    xF = fft(aOut{ii},N,2);
    for nn = 1:N
        sFHat(:,nn)  = hF(:,:,nn) * xF(:,nn);
    end
    sFHat   = sFHat / sqrt(mean(abs(vec(sFHat(:,usedSubcarrierIdx))).^2));
    errMat  = sF{ii}(:,usedSubcarrierIdx) - sFHat(:,usedSubcarrierIdx);
    fVal(ii)= mean(abs(errMat(:)).^2);
end

evmOut = mean(fVal) * 100;