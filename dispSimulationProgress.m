function dispSimulationProgress(opType,arg2)
persistent numIter prevTime
if nargin == 0
    clear numIter prevTime
    return
end

switch lower(opType)
    case 'init'
        dispSimulationProgress();
        numIter = arg2;
        prevTime= 0;
        tic;
    case 'update'
        currIter = arg2;
        currTime = toc;
        if currTime - prevTime > 5
            disp(strcat(['Simulation for this data point is ',num2str(currIter/numIter*100),'% complete. Time remaining for this data point is ', num2str(currTime/currIter * numIter),'s']));
            prevTime = currTime;
        end
end