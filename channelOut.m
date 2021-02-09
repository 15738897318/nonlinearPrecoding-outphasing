function dataOut = channelOut(aMat,hT,usedSubcarrierIdx)

N      = size(aMat,2);
hF     = fft(hT,N,3);
xMat   = 1/sqrt(N) * fft(aMat,N,2);
for nn = 1:numel(usedSubcarrierIdx)
    idx            = usedSubcarrierIdx(nn);
    hFxMat(:,nn)   = hF(:,:,idx) * xMat(:,idx);
end
dataOut= hFxMat;
end