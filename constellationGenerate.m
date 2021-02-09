function [constellation, constellationSize] = constellationGenerate(constellationChoice,normVariance)
% - Generate a vector containing all possible constellation points and the
% modulation object
% - Input -
% constellationChoice - Type of constellation e.g. BPSK, QPSK etc
% normVariance        - (0 or 1) Make constellation unit variance
% - Output -
% constellation       - collection of constellation points
% constellationSize   - number of points in the constellation
% modemObj            - Modem object for the modulate command

if nargin == 0
    error('Insufficient arguments');
end

switch upper(constellationChoice)
    case 'BPSK'
        M           = 2;
    case 'QPSK'
        M           = 4;
    case '16-QAM'
        M           = 16;
    case '64-QAM'
        M           = 64;
    case '256-QAM'
        M           = 256;
    otherwise
        error('Invalid constellation type. Valid constellation types are BPSK,QPSK,8-PAM,16-QAM,64-QAM');
end

constellation     = qammod(0:M-1,M);
constellationSize = length(constellation);

if nargin == 2 && normVariance
    constellation = constellation ./ std(constellation,1);
end

%% - Legacy code %%
% [tempX tempY] = meshgrid(constellationI,constellationQ);
% constellation = tempX +1i* tempY;
% constellation = unique(constellation);