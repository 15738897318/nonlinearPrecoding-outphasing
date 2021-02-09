function [aMat] = nonlinearPrecoding(uF,hT,lambdaAclr,lambdaPen,theoreticalEfficiency,numPsBit,numIter,M,K,N,usedSubcarrierIdx,noiseVar)
% Implementation of Nonlinear algorithm to minimize spectral emissions
% Author: Karthik Upadhya, Stefan Wesemann. Nokia Bell Labs
% Inputs:
% 1. uF - Information symbols in each subcarrier (frequency domain) that
% are to be transmitted.
% 2. hT - Channel impulse response.
% 3. lambdaAclr - Lagrangian weight for ACLR
% 4. lambdaPen - Lagrangian weight for discrete transmit alphabet
% constraint
% 5. theoreticalEfficiency - Theoretical efficiency of outphasing architecture
% 6. numPsBit - Array containing number of phase-shifter bits in each arms
% 7. numIter - Number of iterations
% 8. M - Number of BS antennas
% 9. K - Number of UEs
% 10. N - Number of OFDM subcarriers
% 11. usedSubcarrierIdx - Indices of the non-empty subcarriers.
% 12. noiseVar - Noise variance
unusedSubcarrierIdx = setdiff(1:N,usedSubcarrierIdx);
Nu          = numel(usedSubcarrierIdx);
F           = 1/sqrt(N) * dftmtx(N);
Fz          = F(:,unusedSubcarrierIdx);

% - Initialize to quantized ZF - %
a0          = zfPrecoding(uF,hT,usedSubcarrierIdx,M,N); % Generate zero-forcing precoded signal
aMat        = phaseQuantizeCe(a0,2^max(numPsBit)); % quantize to unit modulus

if theoreticalEfficiency == 1
    alphabet    = exp(1i * 2 * pi * (0:2^max(numPsBit)-1)/2^max(numPsBit)).'; % Constant-envelope precoding with discrete phase-shifts if theoretical efficiency is 1
else
    alphabet    = genAlphabet(numPsBit); % alphabet is non-constant-envelope
end
uT          = sqrt(N) * ifft(uF,N,2); % Generate time-domain OFDM signal

L           = size(hT,3);
hTConcat    = reshape(hT,K,[]);
for nn = 1:N
    if nn >= L
        uTHat(:,nn) = hTConcat * vec(aMat(:,nn:-1:nn-(L-1)));
    else
        uTHat(:,nn) = hTConcat * vec([aMat(:,nn:-1:1) , aMat(:,N:-1:N-(L-nn-1))]);
    end
end
alphaVal = real(uT(:)' * uTHat(:)) / ( norm( uTHat , 'fro' )^2 + K * N * noiseVar ); % Calculate initial value of alpha
rMat     = uT - alphaVal * uTHat; % residual
sqNormAMat = norm(aMat,'fro')^2;

for mm = 1:M
    hmNorm(mm) = norm(vec(squeeze(hT(:,mm,:))))^2;
    hmTildeMat(:,mm) = vec(squeeze(hT(:,mm,:)));
end

aMatFz   = aMat * Fz;

% - Begin non-linear precoding algorithm - %
for ii = 1:numIter
    lambdaPenCurrentIter = lambdaPen^(ii);
    thetaValPremultiplier= 1;
    [~,orderVec]         = sort(abs(aMat(:)),'descend');
    for pp = 1:M*N
        [mm,nn] = ind2sub(size(aMat),orderVec(pp));
        hmTilde = hmTildeMat(:,mm);
        amnBar  = aMat(mm,nn);
        thetaVal= sqrt( max( thetaValPremultiplier * theoreticalEfficiency * M *N - sqNormAMat + abs(amnBar)^2 , 0) );
        fZn     = 1/sqrt(N) * exp(-1i * 2 * pi * (nn-1) * (unusedSubcarrierIdx-1) / N );
        
        if nn > N-L+1
            rBarTilde= vec([rMat(:,nn:end),rMat(:,1:L-(N-nn+1))]);
        else
            rBarTilde= vec(rMat(:,nn:nn+L-1));
        end
        gHg     = alphaVal^2 * hmNorm(mm) + lambdaAclr * (N-Nu)/N + lambdaPenCurrentIter;
        gHY     = alphaVal * hmTilde' * rBarTilde + alphaVal^2 * amnBar * hmNorm(mm) + lambdaAclr * conj(fZn) * (amnBar * fZn.' - aMatFz(mm,:).') + lambdaPenCurrentIter * alphabet.';
        aOpt    = gHY / gHg;
        idx     = abs(aOpt) < thetaVal;
        aOpt(idx) = aOpt(idx) ./ abs(aOpt(idx)) * thetaVal;
        normY   = norm( [rBarTilde + alphaVal * amnBar * hmTilde ; sqrt(lambdaAclr) * (amnBar * fZn.' - aMatFz(mm,:).') ] )^2 + lambdaPenCurrentIter * abs(alphabet.').^2;
        tVal    = normY + abs(aOpt).^2 * gHg - 2 * real( gHY .* conj(aOpt) );
        [minVal,idx]= min(tVal);

        if abs(amnBar) < thetaVal
            fVal = inf;
        else
            fVal = min(normY + abs(amnBar)^2 * gHg - 2 * real( gHY * conj(amnBar) ));
        end
        
        if minVal < fVal % update coordinate only if it reduces the value of the objective function.
            sqNormAMat          = sqNormAMat - abs(amnBar)^2 + abs(aOpt(idx))^2;
            
            if nn > N - L + 1
                rMat(:,nn:end)          = rMat(:,nn:end) - alphaVal * (aOpt(idx) - amnBar) * squeeze(hT(:,mm,1:(N-nn+1)));
                rMat(:,1:L-(N-nn+1))    = rMat(:,1:L-(N-nn+1)) - alphaVal * (aOpt(idx) - amnBar) * squeeze(hT(:,mm,(N-nn+1)+1:end));
            else
                rMat(:,nn:nn+L-1)= rMat(:,nn:nn+L-1) - alphaVal * (aOpt(idx) - amnBar) * squeeze(hT(:,mm,:));
            end
            aMatFz(mm,:)        = aMatFz(mm,:) + (aOpt(idx) - amnBar) * fZn;
            aMat(mm,nn)         = aOpt(idx); % Update coordinate with value that minimizes the objective the most.            
        end
    end
    for nn = 1:N
        if nn >= L
            uTHat(:,nn) = hTConcat * vec(aMat(:,nn:-1:nn-(L-1)));
        else
            uTHat(:,nn) = hTConcat * vec([aMat(:,nn:-1:1) , aMat(:,N:-1:N-(L-nn-1))]);
        end
    end
    alphaVal = real(uT(:)' * uTHat(:)) / ( norm( uTHat , 'fro' )^2 + K * N * noiseVar ); % Calculate initial value of alpha
end
aMat = phaseQuantize(aMat,alphabet); % Quantize the output of the non-linear precoding algorithm to satisfy the discrete alphabet constraints
end


function [output] = phaseQuantize(input,alphabet)
% - Quantize the continuous-input to a given alphabet - %
inputSize   = size(input);
input       = reshape(input,[],1);

dist        = abs(repmat(input,[1,numel(alphabet)]) - repmat(alphabet.',[numel(input),1]));
[~,idx]     = min(dist.');

output      = alphabet(idx);

output     = reshape(output,inputSize);
end

function [output] = phaseQuantizeCe(input,numDacLevel)
% - Quantize the input to a constant-envelope signal with discrete phase
% shifts - %
inputSize   = size(input);
input       = reshape(input,[],1);

alphabet    = exp(1i * 2 * pi * (0:numDacLevel-1)/numDacLevel).';
dist        = abs(repmat(input,[1,numDacLevel]) - repmat(alphabet.',[numel(input),1]));
[~,idx]     = min(dist.');

output      = alphabet(idx);

output      = reshape(output,inputSize);
end

function alphabet = genAlphabet(numPsBit)
% - Generate constellation alphabet at the output of the outphasing
% architecture - %
numDacLevel = 2.^numPsBit;
alphabet1   = exp(1i * 2 * pi * (0:numDacLevel(1)-1).' / numDacLevel(1));
alphabet2   = exp(1i * 2 * pi * (0:numDacLevel(2)-1).' / numDacLevel(2));

alphabet    = 1/2 * vec(alphabet1 + alphabet2.');
cnt = 1;
while(1)
    idx = abs(alphabet(cnt) - alphabet) < 1e-10;
    idx(cnt) = 0;
    alphabet(idx) = [];
    cnt = cnt + 1;
    if cnt > length(alphabet)
        break
    end
end
end