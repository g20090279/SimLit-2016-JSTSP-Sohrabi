function [chnMat,chnAll] =  generate_channel(numAntTx, numAntRx, numPath)
% GENERATE_CHANNEL generates a simplified cluster (also called geometric)
% channel model.
%
% Inputs:
% - numAntTx [scalar]: the number of transmit antennas
% - numAntRx [scalar]: the number of receive antennas
% - numPath [scalar]:  the number of paths
%
% Outputs:
% - chnMat: the generated channel matrix with size numAntRx-by-numAntTx
% - chnAll: all information for this geometric channel model
%

gain = 1/sqrt(2)*randn(numPath,1) + 1i/sqrt(2)*randn(numPath,1); % path gain, following CN(0,1)
phaseTx = 2*pi*rand(numPath,1); % random phase U(0,2pi)
phaseRx = 2*pi*rand(numPath,1);
d = 1/2; % antenna space half wavelength
k = 2*pi;  % 2*pi*lambda*d

% array factor in matrix form (numAnt x numPath)
arrayFactorTx = 1/sqrt(numAntTx)*exp(1i*k*d*(0:numAntTx-1)'*sin(phaseTx.'));
arrayFactorRx = 1/sqrt(numAntRx)*exp(1i*k*d*(0:numAntRx-1)'*sin(phaseRx.'));
g = diag(gain);

% take the number of antennas into account, remove the effect of number of
% paths: the more antennas, the larger power of the channel
chnMat = sqrt(numAntTx*numAntRx/numPath)*arrayFactorRx*g*arrayFactorTx';
chnAll = struct(...
    'pathGain', gain, ...
    'angleDep', phaseTx, ...
    'angleArr', phaseRx, ...
    'antSpace', d, ...
    'numPath', numPath ...
    );

end