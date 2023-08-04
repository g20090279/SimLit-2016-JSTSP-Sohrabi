% Fig. 2: Single-User with Infinite Solution Phase Shifters
% Last revised: August 04, 2023
% Zekai Liang, liang@ice.rwth-aachen.de

%% System Parameters
numAntTx = 64;
numAntRx = 16;
numPath = 15;
numStream = 6;
numRfTx = numStream;
numRfRx = numRfTx;

snrVecDb = (-10:2:6).';
snrVec = 10.^(snrVecDb/10);
snrLen = length(snrVec);

numMC = 100;  % number of iterations for Monte-Carlo simulation

% generate and save channels first for Monte-Carlo simulation
% generate_and_save_channel(numAntTx,numAntRx,numPath,numMC);

% Initial Variables
capFd = 0;     % capacity for fully-digital beamforming
capProp = 0;    % capacity for hybrid beamforming

for iMC = 1:numMC
    %% Generate Channels
    % generate or load one channel instance
    chn = generate_channel(numAntTx,numAntRx,numPath);
    % load(['./data/channels_64x16/channel-',num2str(iMC),'.mat']);
    
    %% Fully-Digital Beamforming
    % SVD Decomposition of the Channel: H = U*D*V'
    gainChn = maxk(eig(chn'*chn),numStream);
    
    % Note 1: Without loss of generality, set the noise power to be 1, so
    % that SNR is equivalent to the power. 
    % Note 2: If the TX and RX beamformers are respectively V and U', the
    % power of the beamformer is equal to the number of streams. Since all
    % the transmit power is defined in the SNR, we need to remove the power
    % of the transmit beamformer, i.e. normalize the transmit beamformer to
    % be unit power.
    capFdTmp = sum( log2(1 + 1/numStream*gainChn*snrVec.') ).';
    capFd = capFd + capFdTmp;

    %% Proposed Algorithm
    capPropTmp = zeros(snrLen,1);
    Vrf = 0;
    for iSnr = 1:snrLen
        % TX Analog Part (Large-Scale Antenna Assumption Vd*Vd~=gamma^2*I)
        F1 = chn' * chn;
        [Vrf,~] = hbf_algorithm1(F1,snrVec(iSnr),numRfTx);
        % Vrf = exp(1i*2*pi*rand(numAntTx,numRfTx)); % DELETE
        
        % TX Digital Part
        Q = Vrf'*Vrf;     % Nrf-by-Nrf
        Heff = chn*Vrf;   % Nt-by-Nrf
        [~,~,V] = svd(Heff*Q^(-1/2));
        
        % Original Model. Assume Equal Power on Each Stream
        Vd = Q^(-1/2)*V(:,1:numStream)/sqrt(numStream);  % normalize to norm(Vrf*Vd,'fro')=1
        
        % A Simplified Model, Assume Equal Power on Each Stream: Vd~=gamma*Ue
        % Vd = V(:,1:numStream)/sqrt(numRfTx*numAntTx);

        % RX Analog Part
        Vt = Vrf*Vd;     % normalized to be 1
        F2 = snrVec(iSnr)*chn*(Vt*Vt')*chn';
        [Wrf,~] = hbf_algorithm1(F2,numAntRx,numRfRx);
        % Wrf = exp(1i*2*pi*rand(numAntRx,numRfRx)); % DELETE
        
        % RX Digital Part, MMSE
        J = snrVec(iSnr)*Wrf'*(chn*(Vt*Vt')*chn')*Wrf + Wrf'*Wrf;
        Wd = sqrt(snrVec(iSnr))*J^(-1)*Wrf'*chn*Vt;

        % The Overall Hybrid Combiner
        Wt = Wrf*Wd; % Note: the power of Wt doesn't affect capacity

        % The Capacity
        capPropTmp(iSnr) = abs(log2(det( eye(numAntRx) + ...
            snrVec(iSnr)*Wt*(Wt'*Wt)^(-1)*Wt'*(chn*(Vt*Vt')*chn') )));
    end

    % Record for Monte-Carlo
    capProp = capProp + capPropTmp;

end

% average to obtain final capacity
capFd = capFd / numMC;
capProp = capProp / numMC;

% plot figure
figure();
grid on; hold on;
plot(snrVecDb,capFd,'k--+',snrVecDb,capProp,'b-o',LineWidth=1);
xlabel('SNR (dB)');
ylabel('Spectral Efficiency (bits/s/Hz)')
lgd = legend('Optimal Fully-Digital Beamforming','Proposed Hybrid Beamforming Algorithm');
lgd.Location = "northwest";
