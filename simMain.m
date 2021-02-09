function simMain(simStart,simEnd)

simData         = readtext('simList.csv');                                  % Read simulation list from simList.csv
varList         = simData(1,:);                                             % Names of variables to be passed to algo
varVals         = simData(2:end,:);                                         % Value of each of the variables
numSim          = size(varVals,1);                                          % Number of simulations to be done
if nargin ~= 2
    simStart        = 1;                                                    % Starting line number of CSV File
    simEnd          = numSim;                                               % Ending line number of CSV File
else
    if ischar(simStart)
        simStart = str2double(simStart);
    end
    if ischar(simEnd)
        simEnd   = str2double(simEnd);
    end
end

if simEnd < simStart || simEnd > numSim
    error('Set simStart and simEnd correctly');
end

for ii = simStart:simEnd
    clc;
    disp(strcat(['************ ***************************************************** *************']));
    disp(strcat(['************ - Running with parameters from row ',num2str(ii),' of simlist.csv - *************']));
    disp(strcat(['************ ***************************************************** *************']));
    if all(cellfun(@isempty,varVals(ii,:)))
        continue; % Skip over empty cells
    elseif strncmpi(varVals{ii,1},'%',1)
        continue; % '%' signifies comment. Skip over commented cells
    end
    
    for jj = 1:length(varList)
        eval([varList{jj} , ' =  varVals{ii,jj};']);
    end
    
    % - Generate directory for dumping results - %
    simSubDir       = sprintf('Results');
    dumpDir         = strcat(fullfile(pwd,simSubDir),filesep);
    if ~isfolder(dumpDir)
        mkdir(dumpDir);
    end
    
    % - Run Algorithm - %
    algo(numTrial,numSubcarrier,numUsedsubcarrier,numBsAntenna,numUe,numChannelTap,numDacBit1,numDacBit2,numIter,lambdaAclr,lambdaPen,theoreticalEfficiency,dumpDir,ii);    
end