function algo(numTrial,numSubcarrier,numUsedsubcarrier,numBsAntenna,numUe,numChannelTap,numDacBit1,numDacBit2,numIter,lambdaAclr,lambdaPen,theoreticalEfficiency,simSubDir,simId)
usedSubcarrierIdx   = (-numUsedsubcarrier/2:numUsedsubcarrier/2-1) + numSubcarrier/2;
constellation       = constellationGenerate('16-QAM',1);

numDacBit           = [numDacBit1,numDacBit2];
noiseVar            = 1;

dispSimulationProgress('init',numTrial);
for ii = 1:numTrial
    dispSimulationProgress('update',ii);
    % - generate channel - %
    hT{ii}                  = complex(normrnd(0,1/sqrt(2),numUe,numBsAntenna,numChannelTap),normrnd(0,1/sqrt(2),numUe,numBsAntenna,numChannelTap));
    
    % - generate payload data - %
    uF{ii}                  = zeros(numUe,numSubcarrier);
    uF{ii}(:,usedSubcarrierIdx) = randsrc(numUe,numUsedsubcarrier,constellation);
    
    % - generate non-linear precoded signal - %
    aMat                    = nonlinearPrecoding(uF{ii},hT{ii},lambdaAclr,lambdaPen,theoreticalEfficiency,numDacBit,numIter,numBsAntenna,numUe,numSubcarrier,usedSubcarrierIdx,noiseVar);
    aOut{ii}                = aMat;
    psdVec(:,ii)            = mean( abs(  fft( aOut{ii} , numSubcarrier , 2 ) ).^2 , 1).'; % calculate PSD
end

save(strcat(simSubDir,num2str(simId),'.mat'),'hT','aOut','uF');