% Multiuser
% Last revised: Mar 15, 2023
% Zekai Liang, liang@ice.rwth-aachen.de

%% Settings
snrVecDb = (-10:10).';
snrVec = 10.^(snrVecDb/10);
snrLen = length(snrVec);

addpath('C:\Users\liang\Documents\Git-Local\gitlab-iss\communication-systems');
PA = PowerAllocator();

%% RF precoder
numUser = 8;
numAntTx = 64;
numAntRx = 1;
numRfTx = numUser + 1;
numRfTx2 = numUser;
numPath = 15;

%% Generate or Load Channels
numMC = 100;  % number of iterations for Monte-Carlo simulation
% generate_and_save_channel(numAntTx,numAntRx,numPath,numMC);\

%% Initialize Variables for Monte-Carlo Simulation
capFdSum = 0;
capPropSum = 0;

for iMC = 1:numMC

    %% Load or Generate Channels
    % Generate Multiuser Channel Instance
%     generate_and_save_channel_mu(numUser,numAntTx,numAntRx,numPath,numMC,'combined');

    % Load Channel Data
    load(['./data/channels_mu_64x1/channel-',num2str(iMC),'.mat']);
    
    %% Fully-Digital Precoder
    Vd = chn'*(chn*chn')^(-1);
    
    % For Each SNR
    capFd = zeros(numUser,snrLen);
    for iSnr = 1:snrLen
        % Use Water-Filling to Allocate Power
        powFd = PA.water_filling_MIMO_ZF(snrVec(iSnr),ones(numUser,1),Vd);
        VdNew = Vd*diag(sqrt(powFd));
    
        % Calculate Multi-User Capacity, i.e. Sum-Rate of All Users
        for iUser = 1:numUser
            chnUser = chn(iUser,:);
            gainTmp = abs(chnUser*VdNew).^2;
            gainTmpBar = gainTmp;
            gainTmpBar(iUser) = [];
            capFd(iUser,iSnr) = log2( 1 + gainTmp(iUser) / (1 + sum(gainTmpBar)) );
        end
    end
    
    % the Sum Capacity of All Users
    capFdSum = capFdSum + sum(capFd)';


    %% Proposed Algorithm

    % Caculate Analog Precoder by Algorithm 3
    Vrf = hbf_algorithm3(chn,numAntTx,numRfTx);
    
    % Caculate Digital Precoder with Zero-Forcing
    Vd = Vrf'*chn'*(chn*(Vrf*Vrf')*chn')^(-1);

    % For Each SNR
    capMu = zeros(numUser,snrLen);
    
    for iSnr = 1:snrLen
       
        % Water-Filling Algorithm for Power Allocation
        % Note that for low-SNR case, not every stream is allocated power.
        powProp = PA.water_filling_MIMO_ZF(snrVec(iSnr),ones(numUser,1),Vrf*Vd);
        VdNew = Vd*sqrt(diag(powProp));
        
        % Calculate Multi-User Capacity, i.e. Sum-Rate of All Users
        for iUser = 1:numUser
            chnUser = chn(iUser,:);
            gainTmp = abs(chnUser*Vrf*VdNew).^2;
            gainTmpBar = gainTmp;
            gainTmpBar(iUser) = [];
            capMu(iUser,iSnr) = log2( 1 + gainTmp(iUser) / (1 + sum(gainTmpBar)) );
        end

    end

    % the Sum Capacity of All Users
    capPropSum = capPropSum + sum(capMu)';

    %
    disp([num2str(iMC),'-th iteration finished!'])
end

%% Averge Data After Monte-Carlo
capFdSum = capFdSum/numMC;
capPropSum = capPropSum/numMC;

%% Plot Figure
figure;
plot(snrVecDb,capFdSum,'--k',snrVecDb,capPropSum,'-ob');
hold on;
grid on;
xlabel('SNR (dB)');
ylabel('Sum Rate (bits/s/Hz)');
legend("Fully-Digital ZF","Proposed Algorithm for N^{RF}=9");
title('Figure 5')