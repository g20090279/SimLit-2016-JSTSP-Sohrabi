% Multiuser
% Last revised: June 09, 2023
% Zekai Liang, liang@ice.rwth-aachen.de

%% Settings
snrVecDb = (-10:10).'; % set high SNR to check water-filling algorithm and algorithm 3
snrVec = 10.^(snrVecDb/10);
snrLen = length(snrVec);

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
capRef33Sum = 0;
capRef32_1Sum = 0;
capRef32_2Sum = 0;
capRef32Sum = 0;

for iMC = 1:numMC

    %% Load or Generate Channels

    % NOTE: the dimension of generated channel of each user is a row
    % vector. Therefore, the combined channel matrix is H=[h1',...,hK']',
    % where (.)' is the conjugate transpose (Hermitian).

    % Generate Multiuser Channel Instance
%     generate_and_save_channel_mu(numUser,numAntTx,numAntRx,numPath,numMC,'combined');

    % Load Channel Data
    load(['./data/channels_mu_64x1/channel-',num2str(iMC),'.mat']);
    
    %% Fully-Digital Precoder
    Vd = chnMat'*(chnMat*chnMat')^(-1);
    
    % For Each SNR (Assume the Noise Power equals 1)
    capFd = zeros(numUser,snrLen);
    for iSnr = 1:snrLen

        % Use Water-Filling to Allocate Power, to Meet the Power Constraint
        % (The total beamformer power constraint is ||Vd||_F^2=P)
        powFd = water_filling_MIMO_ZF(snrVec(iSnr),ones(numUser,1),Vd);
        VdNew = Vd*diag(sqrt(powFd));
    
        % Calculate Multi-User Capacity, i.e. Sum-Rate of All Users
        for iUser = 1:numUser
            chnUser = chnMat(iUser,:);
            powerRx = abs(chnUser*VdNew).^2;
            capFd(iUser,iSnr) = log2( 1 + powerRx(iUser) / ...
                ( 1 + sum(powerRx([1:iUser-1,iUser+1:numUser])) ) );
        end
    end
    
    % the Sum Capacity of All Users
    capFdSum = capFdSum + sum(capFd)';


    %% Proposed Algorithm: Iterative Analog + ZF Digital
    % For Each SNR (Assume the Noise Power equals 1)
    capMu = zeros(numUser,snrLen);
    for iSnr = 1:snrLen

        % Caculate Analog Precoder by Algorithm 3
        [Vrf,Vd,powProp] = hbf_algorithm3(chnMat,numAntTx,numRfTx,numUser,snrVec(iSnr));
        
        % Add Power Factor to Digital Precoder
        VdNew = Vd*sqrt(diag(powProp));
        
        % Calculate Multi-User Capacity, i.e. Sum-Rate of All Users
        for iUser = 1:numUser
            chnUser = chnMat(iUser,:);
            powerRx = abs(chnUser*Vrf*VdNew).^2;
            capMu(iUser,iSnr) = log2( 1 + powerRx(iUser) / ...
                ( 1 + sum(powerRx([1:iUser-1,iUser+1:numUser])) ) );
        end

    end

    % the Sum Capacity of All Users
    capPropSum = capPropSum + sum(capMu)';

    
    %% Reference [33]: Analog Conjugate Transpose (eqv. MRC) + Digital ZF
    % The Phase of Analog Beamformer Is the Phase of Conjugate Transpose of
    % the Composite Channel Matrix. This Is Equivalent to MRC. But In MU
    % Case, After MRC Analog Beamformer, Interuser-Interference (IUI) Still
    % Exists. The Digital ZF Helps to Eliminate the IUI.
    Vrf = 1/numAntTx * exp(1i*angle(chnMat'));
    
    % the Equivalent Channel
    chnEqv = chnMat*Vrf;

    % ZF Digital Beamformer
    Vd = chnEqv'*(chnEqv*chnEqv')^-1;

    % For Each SNR (Assume the Noise Power equals 1)
    capRef33 = zeros(numUser,snrLen);
    for iSnr = 1:snrLen
        % Use Water-Filling Algorithm to Meet Power Constraint
        % (The total beamformer power constraint is ||Vrf*Vd||_F^2=P)
        powAlloc = water_filling_MIMO_ZF(snrVec(iSnr),ones(numUser,1),Vrf*Vd);
        VdNew = Vd*diag(sqrt(powAlloc));

        % Calculate Multi-User Capacity, i.e. Sum-Rate of All Users
        for iUser = 1:numUser
            chnUser = chnMat(iUser,:);
            powerRx = abs(chnUser*Vrf*VdNew).^2;
            capRef33(iUser,iSnr) = log2( 1 + powerRx(iUser) / ...
                ( 1 + sum(powerRx([1:iUser-1,iUser+1:numUser])) ) );
        end
    end

    % the Sum Capacity of All Users
    capRef33Sum = capRef33Sum + sum(capRef33)';
    

    %% Reference [32]: Hybrid Ananlog Point to Best Path (MUBS) + Digital ZF
    % Analog Beamformer Pointing to Best Path
    Vrf = zeros(numAntTx, numUser);
    for iUser = 1:numUser
        [~, indMaxGain] = max(abs(chnAll(iUser).pathGain));
        Vrf(:,iUser) = 1/sqrt(numAntTx)*exp(1i*2*pi*1/2*(0:numAntTx-1)'*sin(chnAll(iUser).angleDep(indMaxGain)));
    end

    % the Equivalent Channel
    chnEqv = chnMat*Vrf;

    % ZF Digital Beamformer
    Vd = chnEqv'*(chnEqv*chnEqv')^-1;

    % For Each SNR (Assume the Noise Power equals 1)
    capRef32 = zeros(numUser,snrLen);
    for iSnr = 1:snrLen
        % Use Water-Filling Algorithm to Meet Power Constraint
        % (The total beamformer power constraint is ||Vrf*Vd||_F^2=P)
        powAlloc = water_filling_MIMO_ZF(snrVec(iSnr),ones(numUser,1),Vrf*Vd);
        VdNew = Vd*diag(sqrt(powAlloc));

        % Calculate Multi-User Capacity, i.e. Sum-Rate of All Users
        for iUser = 1:numUser
            chnUser = chnMat(iUser,:);
            powerRx = abs(chnUser*Vrf*VdNew).^2;
            capRef32(iUser,iSnr) = log2( 1 + powerRx(iUser) / ...
                ( 1 + sum(powerRx([1:iUser-1,iUser+1:numUser])) ) );
        end
    end

    % the Sum Capacity of All Users
    capRef32Sum = capRef32Sum + sum(capRef32)';
    

    %% Print the Progress
    disp([num2str(iMC),'-th Monte-Carlo iteration finished!'])

end

%% Averge Data After Monte-Carlo
capFdSum = capFdSum/numMC;
capPropSum = capPropSum/numMC;
capRef33Sum = capRef33Sum/numMC;
capRef32Sum = capRef32Sum/numMC;

%% Plot Figure
figure;
plot(snrVecDb,capFdSum,'--k',snrVecDb,capPropSum,'-ob',...
    snrVecDb,capRef33Sum,'-^r',snrVecDb,capRef32Sum,'-vg');
hold on;
grid on;
xlabel('SNR (dB)');
ylabel('Sum Rate (bits/s/Hz)');
legend("Fully-Digital ZF","Proposed Algorithm for N^{RF}=9",...
    "Hybrid beamforming in [33], N^{RF}=8", ...
    "Hybrid Beamforming in [32], N^{RF}=8");
