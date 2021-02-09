clc
clear all
close all
simData         = readtext('simList.csv');                                  % Read simulation list from simList.csv
varList         = simData(1,:);                                             % Names of variables to be passed to algo
varVals         = simData(2:end,:);                                         % Value of each of the variables
numSim          = size(varVals,1);                                          % Number of simulations to be done

% - Plot 1 - %
xData = [];
yData = [];
for ii = 1:numSim
    simID       = ii;
    if all(cellfun(@isempty,varVals(ii,:)))
        continue; % Skip over empty cells
    elseif strncmpi(varVals{ii,1},'%',1)
        continue; % '%' signifies comment. Skip over commented cells
    end
    
    for jj = 1:length(varList)
        eval([varList{jj} , ' =  varVals{ii,jj};']);
    end
    
    if theoreticalEfficiency ~= 1 || numBsAntenna ~= 32
        continue
    end
    try
        load(strcat('Results/',num2str(ii),'.mat'));
    catch
        continue
    end
    yData = [yData aclrCalc(aOut,256)];
    xData = [xData evmCalc(uF,hT,aOut,256)];
end
figure(1); plot(xData,yData,'b*-','DisplayName','Constant-envelope')
hold all
[xData,psdData] = psdPlot(aOut,numTrial,numSubcarrier,numUsedsubcarrier);
figure(2); plot(xData,psdData,'b-','DisplayName','Constant-envelope');
hold all

xData = [];
yData = [];
for ii = 1:numSim
    simID       = ii;
    if all(cellfun(@isempty,varVals(ii,:)))
        continue; % Skip over empty cells
    elseif strncmpi(varVals{ii,1},'%',1)
        continue; % '%' signifies comment. Skip over commented cells
    end
    
    for jj = 1:length(varList)
        eval([varList{jj} , ' =  varVals{ii,jj};']);
    end
    
    if theoreticalEfficiency ~= 0.75 || numBsAntenna ~= 32
        continue
    end
    try
        load(strcat('Results/',num2str(ii),'.mat'));
    catch
        continue
    end
    yData = [yData aclrCalc(aOut,256)];
    xData = [xData evmCalc(uF,hT,aOut,256)];
end
figure(1); plot(xData,yData,'ro-','DisplayName','\eta = 0.75');

[xData,psdData] = psdPlot(aOut,numTrial,numSubcarrier,numUsedsubcarrier);
figure(2); plot(xData,psdData,'r-','DisplayName','\eta = 0.75');


xData = [];
yData = [];
for ii = 1:numSim
    simID       = ii;
    if all(cellfun(@isempty,varVals(ii,:)))
        continue; % Skip over empty cells
    elseif strncmpi(varVals{ii,1},'%',1)
        continue; % '%' signifies comment. Skip over commented cells
    end
    
    for jj = 1:length(varList)
        eval([varList{jj} , ' =  varVals{ii,jj};']);
    end
    
    if theoreticalEfficiency ~= 0.5 || numBsAntenna ~= 32
        continue
    end
    
    try
        load(strcat('Results/',num2str(ii),'.mat'));
    catch
        continue
    end
    yData = [yData aclrCalc(aOut,256)];
    xData = [xData evmCalc(uF,hT,aOut,256)];
end
figure(1); plot(xData,yData,'k^-','DisplayName','\eta = 0.5')

[xData,psdData] = psdPlot(aOut,numTrial,numSubcarrier,numUsedsubcarrier);
figure(2); plot(xData,psdData,'k-','DisplayName','\eta = 0.5');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xData = [];
yData = [];
for ii = 1:numSim
    simID       = ii;
    if all(cellfun(@isempty,varVals(ii,:)))
        continue; % Skip over empty cells
    elseif strncmpi(varVals{ii,1},'%',1)
        continue; % '%' signifies comment. Skip over commented cells
    end
    
    for jj = 1:length(varList)
        eval([varList{jj} , ' =  varVals{ii,jj};']);
    end
    
    if theoreticalEfficiency ~= 1 || numBsAntenna ~= 128
        continue
    end
    try
        load(strcat('Results/',num2str(ii),'.mat'));
    catch
        continue
    end
    yData = [yData aclrCalc(aOut,256)];
    xData = [xData evmCalc(uF,hT,aOut,256)];
end
figure(1); plot(xData,yData,'b*--','DisplayName','Constant-envelope')
hold all

xData = [];
yData = [];
for ii = 1:numSim
    simID       = ii;
    if all(cellfun(@isempty,varVals(ii,:)))
        continue; % Skip over empty cells
    elseif strncmpi(varVals{ii,1},'%',1)
        continue; % '%' signifies comment. Skip over commented cells
    end
    
    for jj = 1:length(varList)
        eval([varList{jj} , ' =  varVals{ii,jj};']);
    end
    
    if theoreticalEfficiency ~= 0.75 || numBsAntenna ~= 128
        continue
    end
    try
        load(strcat('Results/',num2str(ii),'.mat'));
    catch
        continue
    end
    yData = [yData aclrCalc(aOut,256)];
    xData = [xData evmCalc(uF,hT,aOut,256)];
end
figure(1); plot(xData,yData,'ro--','DisplayName','\eta = 0.75')
xData = [];
yData = [];
for ii = 1:numSim
    simID       = ii;
    if all(cellfun(@isempty,varVals(ii,:)))
        continue; % Skip over empty cells
    elseif strncmpi(varVals{ii,1},'%',1)
        continue; % '%' signifies comment. Skip over commented cells
    end
    
    for jj = 1:length(varList)
        eval([varList{jj} , ' =  varVals{ii,jj};']);
    end
    
    if theoreticalEfficiency ~= 0.5 || numBsAntenna ~= 128
        continue
    end
    
    try
        load(strcat('Results/',num2str(ii),'.mat'));
    catch
        continue
    end
    yData = [yData aclrCalc(aOut,256)];
    xData = [xData evmCalc(uF,hT,aOut,256)];
end
figure(1); plot(xData,yData,'k^--','DisplayName','\eta = 0.5')

legend('show')
xlabel('EVM (%)');
ylabel('ACLR (dB)');

figure(2); 
legend('show')
xlabel('Subcarrier index');
ylabel('ACLR (dB)');
