function [BF,cap] = hbf_algorithm1(F,snr,numAnt,numRf,numStream)
% Convergence threshold
epsilon = 10^-3;

% Initialize analog beamformer with size #ant-by-#RFchain
BF = 1/sqrt(numAnt)*exp(1i*2*pi*rand([numAnt,numRf]));

% Calculate initial capacity
% Note that it doesn't matter if we normalize analog beamformer with
% 1/sqrt(numRf) or not. Since later we have to obtain digital beamformer,
% and there we will assure the combined hybrid beamformer has unit power.
cap = abs(log2(det( eye(numRf)+snr*BF'*F*BF )));

% Obtain transmit analog beamformer

% Reset variables
diff = inf;
count = 0;

while diff > epsilon
    for rfInd = 1:numRf
        VrfBar = BF;
        VrfBar(:,rfInd) = [];
        Cj = eye(numStream-1) + snr*VrfBar'*F*VrfBar;
        Gj = snr*F-snr^2*F*VrfBar*Cj^(-1)*VrfBar'*F;
        for aInd = 1:numAnt
            gj = Gj(aInd,:);
            gj(rfInd) = [];
            vrfj = BF(:,rfInd);
            vrfj(rfInd) = [];
            BF(aInd,rfInd) = 1/sqrt(numAnt)*exp(1i*angle(gj*vrfj));  % note: update the current Vrf
        end
    end

    % Calculate new capacity, assuming Vd*Vd'=I
    capTxTmp = abs(log2(det( eye(numRf) + snr*BF'*F*BF )));
    diff = abs(capTxTmp-cap)/capTxTmp;
    cap = capTxTmp;

    % count how many iterations to converge
    count = count+1;
end

% record the final capacity
disp(['SNR=',num2str(10*log10(snr)),'dB: converged after ',num2str(count),' iterations.']);
end