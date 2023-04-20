% Last revised: Mar 15, 2023
% Zekai Liang, liang@ice.rwth-aachen.de

%% Fig 2: SU, capacity, infinite solution
numAntTx = 64;
numAntRx = 16;
numPath = 15;
numStream = 6;
numRfTx = numStream; % recommend to design Nrf=Ns even Ns<Nrf<2Ns
numRfRx = numRfTx;

snrVecDb = (-10:6).';
snrVec = 10.^(snrVecDb/10);
snrLen = length(snrVec);

numMC = 100;  % number of iterations for Monte-Carlo simulation

% generate and save channels first for Monte-Carlo simulation
% generate_and_save_channel(numAntTx,numAntRx,numPath,numMC);

% initial variables
capSVD = 0;     % capacity for fully-digital beamformer
capPropTx = 0;  % TX analog
capProp = 0;    % TX + RX
capPropRx = 0;  % capacity only RX beamforanalog part

for iMC = 1:numMC
    % Generate Channels
    % generate or load one channel instance
    % chn = generate_channel(numAntTx,numAntRx,numPath)
    load(['./data/channels/channel-',num2str(iMC),'.mat']);
    
    %%% Fully-Digital Beamforming, No Water-Filling
    [~,svVal,V] = svd(chn);  % transmit fully-digital beamformer V*V'~=I!
    svVal = diag(svVal);
    svdGain = svVal(1:numStream).^2;
    
    % check the pwer of the beamformer
    bfDigital = V(:,1:numStream);
    powBfDigital = trace(bfDigital'*bfDigital);  % equals to numStream, since V is unitary
    
    % Since we consider the SNR, we have to make P=1. Therefore, we need to
    % compensate with a power factor to make the power of beamformer to be 1,
    % i.e., add factor 1/numStream and make sure Tr(svdGain*svdGain')=1
    capSvdThis = log2(1 + 1/numStream*svdGain*snrVec.');
    capSvdThis = sum(capSvdThis).';
    capSVD = capSVD + capSvdThis;

    %%% Proposed Algorithm
    capTxTmp = zeros(snrLen,1);
    capRxTmp = zeros(snrLen,1);
    capPropTmp = zeros(snrLen,1);

    for snrInd = 1:snrLen
        % TX analog part by assuming Vd*Vd~=gamma^2*.I
        F1 = chn' * chn;
        [Vrf,capTxTmp(snrInd)] = hbf_algorithm1(F1,snrVec(snrInd),numAntTx,numRfTx,numStream);
        
        % TX, no water-filling (assume equal power)
        Q = Vrf'*Vrf;     % Nrf-by-Nrf
        Heff = chn*Vrf;   % Nt-by-Nrf
        [~,~,V] = svd(Heff*Q^(-1/2));
        Vd = Q^(-1/2)*V(:,1:numStream)/sqrt(numStream);  % normalize to trace(Vd*Vrf'*Vrf*Vd)=1

        % RX analog part with hair in the toilent
        Vt = Vrf*Vd;
        F2 = chn*(Vt*Vt')*chn';
        [Wrf,capRxTmp(snrInd)] = hbf_algorithm1(F2,snrVec(snrInd),numAntRx,numRfRx,numStream);
        
        % RX digital part, MMSE
        J = Wrf'*(chn*(Vt*Vt')*chn')*Wrf+1/snrVec(snrInd)*(Wrf'*Wrf);
        Wd = J^(-1)*Wrf'*chn*Vt;

        % Calculate the capacity
        Wt = Wrf*Wd;
%         Wt = Wt/norm(Wt,'fro')*sqrt()
        capPropTmp(snrInd) = abs(log2(det( eye(numAntRx)+snrVec(snrInd)*Wt*(Wt'*Wt)^(-1)*Wt'*(chn*(Vt*Vt')*chn') )));
    end

    % record for Monte-Carlo
    capPropTx = capPropTx + capTxTmp;
    capPropRx = capPropRx + capRxTmp;
    capProp = capProp + capPropTmp;

end

% average to obtain final capacity
capSVD = capSVD / numMC;
capProp = capProp / numMC;

% plot figure
figure();
grid on; hold on;
plot(snrVecDb,capSVD,'k--',snrVecDb,capProp,'b-o');
xlabel('SNR (dB)');
ylabel('Spectral Efficiency (bits/s/Hz)')
legend('Optimal Fully-Digital Beamforming','Proposed Hybrid Beamforming Algorithm')


%% MU
