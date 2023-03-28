% Last revised: Mar 15, 2023
% Zekai Liang, liang@ice.rwth-aachen.de

%%% Fig 2: capacity, infinite solution
numAntTx = 64;
numAntRx = 16;
numPath = 15;
numStream = 6;
numRfTx = numStream; % recommend to design Nrf=Ns even Ns<Nrf<2Ns

chn = generate_channel(numAntTx,numAntRx,numPath);

snrVecDb = (-10:6);
snrVec = 10.^(snrVecDb/10);
snrLen = length(snrVec);

%% fully-digital beamforming, no water-filling
[~,svVal,V] = svd(chn);  % transmit fully-digital beamformer V*V'~=I!
svVal = diag(svVal);
svdGain = svVal(1:numStream).^2;

% check the pwer of the beamformer
bfDigital = V(:,1:numStream);
powBfDigital = trace(bfDigital'*bfDigital);  % equals to numStream, since V is unitary

% Since we consider the SNR, we have to make P=1. Therefore, we need to
% compensate with a power factor to make the power of beamformer to be 1,
% i.e., add factor 1/numStream and make sure Tr(svdGain*svdGain')=1
capSVD = log2(1 + 1/numStream*svdGain*snrVec);
capSVD = sum(capSVD).';

%% proposed
F1 = chn'*chn;
epsilon = 10^-3;

% initialize transmit analog beamformer
Vrf = exp(1i*2*pi*rand([numAntTx,numRfTx]));

% obtain transmit digital beamformer, no water-filling (assume equal power)
Q = Vrf'*Vrf;     % Nrf-by-Nrf
Heff = chn*Vrf;   % Nt-by-Nrf
[U,D,V] = svd(Heff*Q^(-1/2));
Vd = Q^(-1/2)*V(:,1:numStream)/sqrt(numStream); % normalize!!!

% calculate initial capacity
capProp = zeros(snrLen,1);
capPropNew = zeros(snrLen,1);
for snrInd = 1:snrLen
    capPropNew(snrInd) = abs(log2(det( eye(numAntRx)+snrVec(snrInd)*Heff*(Vd*Vd')*Heff' )));
end

% obtain transmit analog beamformer
for snrInd = 1:snrLen
    VrfNew = zeros(numAntTx,numRfTx);
    diff = inf;
    count = 0;
    while diff > epsilon
        capProp(snrInd) = capPropNew(snrInd);
        for rfInd = 1:numRfTx
            VrfBar = Vrf;
            VrfBar(:,rfInd) = [];
            Cj = eye(numStream-1) + snrVec(snrInd)*VrfBar'*F1*VrfBar;
            Gj = snrVec(snrInd)*F1-(snrVec(snrInd))^2*F1*VrfBar*Cj^(-1)*VrfBar'*F1;
            for tInd = 1:numAntTx
                gj = Gj(tInd,:);
                gj(rfInd) = [];
                vrfj = Vrf(:,rfInd);
                vrfj(rfInd) = [];
                VrfNew(tInd,rfInd) = exp(1i*angle(gj*vrfj));
            end
        end
        % calculate new digital beamformer wrt the new analog beamformer
        Q = VrfNew'*VrfNew;
        Heff = chn*VrfNew;   % Nt-by-Nrf
        [U,D,V] = svd(Heff*Q^(-1/2));
        Vd = Q^(-1/2)*V(:,1:numStream)/sqrt(numStream);

        % calculate new capacity
        capPropNew(snrInd) = abs(log2(det( eye(numAntRx)+snrVec(snrInd)*Heff*(Vd*Vd')*Heff' )));
        diff = (capPropNew(snrInd)-capProp(snrInd))/capPropNew(snrInd);
    
        % count how many iterations to converge
        count = count+1;
    end

    % record the final capacity
    capProp(snrInd) = capPropNew(snrInd);
    disp(['SNR-',num2str(10*log10(snrVec(snrInd))),'dB converged after ',num2str(count),' iterations.']);
end

% plot figure
figure();
grid on;
plot(snrVecDb,capSVD,'k--',snrVecDb,capProp,'b-o');
xlabel('SNR (dB)');
ylabel('Spectral Efficiency (bits/s/Hz)')
