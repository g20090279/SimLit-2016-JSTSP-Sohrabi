% Generate channel
% Last revised: Mar. 23, 2023
% Zekai Liang, liang@ice.rwth-aachen.de
function channel =  generate_channel(numAntTx, numAntRx, numPath)
% Inputs:
% - numAntTx: the number of transmit antennas
% - numAntRx: the number of receive antennas
% - numPath:  the number of paths
% Outputs:
% - channel: the channel matrix

gain = 1/sqrt(2)*randn(numPath,1) + 1i/sqrt(2)*randn(numPath,1);
phaseTx = 2*pi*rand(numPath,1);
phaseRx = 2*pi*rand(numPath,1);
d = 1/2;
k = 2*pi;  % 2*pi*lambda*d

% array factor in matrix form (numAnt x numPath)
arrayFactorTx = 1/sqrt(numAntTx)*exp(1i*k*d*(0:numAntTx-1)'*sin(phaseTx.'));
arrayFactorRx = 1/sqrt(numAntRx)*exp(1i*k*d*(0:numAntRx-1)'*sin(phaseRx.'));
g = diag(gain);

% take the number of antennas into account, remove the effect of number of
% paths: the more antennas, the larger power of the channel
channel = sqrt(numAntTx*numAntRx/numPath)*arrayFactorRx*g*arrayFactorTx';

end