function a0 = zfPrecoding(sF,hT,usedSubcarrierIdx,M,N)

hF  = fft(hT,N,3);
Nu  = numel(usedSubcarrierIdx); % Number of non-empty subcarriers
% - generate precoded symbols - %
u0F = zeros(M,N);
for nn = 1:Nu
    subcarrierIdx = usedSubcarrierIdx(nn);
    u0F(:,subcarrierIdx) = hF(:,:,subcarrierIdx)' / (hF(:,:,subcarrierIdx) * hF(:,:,subcarrierIdx)') * sF(:,subcarrierIdx); % ZF precoding
end
a0  = sqrt(N) * ifft(u0F,N,2); % Generate precoded symbols in time domain