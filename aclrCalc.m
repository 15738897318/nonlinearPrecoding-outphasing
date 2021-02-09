function aclrOut = aclrCalc(aOut,numUsedsubcarrier,usedSubcarrierIdx)
% - Calculate ACLR at the output of the precoder - %
numTrial = numel(aOut);
N        = size(aOut{1},2);

if ~exist('usedSubcarrierIdx','var')
    usedSubcarrierIdx   = (-numUsedsubcarrier/2:numUsedsubcarrier/2-1) + N/2;
end
unusedSubcarrierIdx     = setdiff(1:N,usedSubcarrierIdx);
    
for ii = 1:numTrial
    psdVec(:,ii) = mean(abs(fft(aOut{ii},N,2)).^2,1);
end
psdOut = mean(psdVec,2);
psdOut = psdOut / mean(psdOut(usedSubcarrierIdx));

aclrOut= 10 * log10( mean(psdOut(unusedSubcarrierIdx)) );