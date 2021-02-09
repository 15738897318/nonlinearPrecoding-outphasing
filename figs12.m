function [xData,psdOut] = psdPlot(aOut,numTrial,numSubcarrier,numUsedsubcarrier)

usedSubcarrierIdx   = (-numUsedsubcarrier/2:numUsedsubcarrier/2-1) + numSubcarrier/2;
   
for ii = 1:numTrial
    psdVec(:,ii) = mean(abs(fft(aOut{ii},numSubcarrier,2)).^2,1);
end
psdOut = mean(psdVec,2);
psdOut = 10*log10(psdOut / mean(psdOut(usedSubcarrierIdx)));

xData = (-(numSubcarrier-1)/2:(numSubcarrier-1)/2);