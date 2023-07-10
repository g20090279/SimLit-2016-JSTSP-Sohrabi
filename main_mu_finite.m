% Multiuser with Finite-Resolution Phase Shifters
% Last revised: July 08, 2023
% Zekai Liang, liang@ice.rwth-aachen.de

%% Settings
% Considered SNR Range
snrVecDb = (-10:10).'; % set high SNR to check water-filling algorithm and algorithm 3
snrVec = 10.^(snrVecDb/10);
snrLen = length(snrVec);

% System Setting
numUser = 4;
numAntTx = 64;
numAntRx = 1;
numRfTx = numUser + 1; % number of RF chain in the proposed algorithm
numPath = 15;
numBit = 1; % considers only very low resolution phase shifter
codeBook = exp(1i*2*pi*(0:1)/(2^numBit)).';

% Simulation Setting
numMC = 100;  % number of iterations for Monte-Carlo simulation

%% Initialize Variables for Monte-Carlo Simulation
capFdSum = 0;
capPropSum = 0;
capPropQuan1Sum = 0;  % quantize final results from infinite resolution
capPropQuan2Sum = 0;  % quantize intermediate results at every iteration
capRef33Sum = 0;
capRef33QuanSum = 0;
capRef32Sum = 0;
capRef32QuanSum = 0;

for iMC = 1:numMC

    %% Load or Generate Channels

    % NOTE: the dimension of generated channel of each user is Nr-by-Nt.
    % Therefore, the MISO channel for each user is a row-vector (h') and
    % the combined channel matrix is H=[h1',...,hK']', where (.)' is the
    % conjugate transpose (Hermitian).

    % Generate Multiuser Channel Instance
%     generate_and_save_channel_mu(numUser,numAntTx,numAntRx,numPath,numMC,'combined');

    % Load Channel Data
    load(['./data/channels_mu_64x1/channel-',num2str(iMC),'.mat']);

    %% User Previously Generated Channel Matrix
    chnMat = chnMat(1:numUser,:);

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
    capProp = zeros(numUser,snrLen);
    capPropQuan1 = zeros(numUser,snrLen);
    capPropQuan2 = zeros(numUser,snrLen);
    
    for iSnr = 1:snrLen

        % Caculate Analog Precoder by Algorithm 3
        [Vrf,Vd,powProp] = hbf_algorithm3(chnMat,numAntTx,numRfTx,...
            numUser,snrVec(iSnr),'infinite');
        [VrfQuan2,VdQuan2,powPropQuan2] = hbf_algorithm3(chnMat,numAntTx,...
            numRfTx,numUser,snrVec(iSnr),'finite',codeBook);

        % Quantize 1: the Analog Precoder to the Nearest Points in Codebook
        VrfQuan1 = quantizeByCodebook(codeBook,Vrf);
        VdQuan1 = VrfQuan1'*chnMat'*(chnMat*(VrfQuan1*VrfQuan1')*chnMat')^(-1);
        powPropQuan1 = water_filling_MIMO_ZF(snrVec(iSnr),ones(numUser,1),VrfQuan1*VdQuan1);
        
        % Add Power Factor to Digital Precoder
        VdNew = Vd*sqrt(diag(powProp));
        VdQuan1New = VdQuan1*sqrt(diag(powPropQuan1));
        VdQuan2New = VdQuan2*sqrt(diag(powPropQuan2));
        
        % Calculate Multi-User Capacity, i.e. Sum-Rate of All Users
        for iUser = 1:numUser
            chnUser = chnMat(iUser,:);

            powerRx = abs(chnUser*Vrf*VdNew).^2;
            powerRxQuan1 = abs(chnUser*VrfQuan1*VdQuan1New).^2;
            powerRxQuan2 = abs(chnUser*VrfQuan2*VdQuan2New).^2;
            
            capProp(iUser,iSnr) = log2( 1 + powerRx(iUser) / ...
                ( 1 + sum(powerRx([1:iUser-1,iUser+1:numUser])) ) );
            capPropQuan1(iUser,iSnr) = log2( 1 + powerRxQuan1(iUser) / ...
                ( 1 + sum(powerRxQuan1([1:iUser-1,iUser+1:numUser])) ) );
            capPropQuan2(iUser,iSnr) = log2( 1 + powerRxQuan2(iUser) / ...
                ( 1 + sum(powerRxQuan2([1:iUser-1,iUser+1:numUser])) ) );
        end

    end

    % the Sum Capacity of All Users
    capPropSum = capPropSum + sum(capProp)';
    capPropQuan1Sum = capPropQuan1Sum + sum(capPropQuan1)';
    capPropQuan2Sum = capPropQuan2Sum + sum(capPropQuan2)';

    %% Reference [33]: Analog Conjugate Transpose (eqv. MRC) + Digital ZF
    % The Phase of Analog Beamformer Is the Phase of Conjugate Transpose of
    % the Composite Channel Matrix. This Is Equivalent to MRC. But In MU
    % Case, After MRC Analog Beamformer, Interuser-Interference (IUI) Still
    % Exists. The Digital ZF Helps to Eliminate the IUI.
    Vrf = 1/numAntTx * exp(1i*angle(chnMat'));

    % Quantize the Analog Beamformer
    VrfQuan1 = quantizeByCodebook(codeBook,Vrf);
    
    % the Equivalent Channel
    chnEqv = chnMat*Vrf;
    chnEqvQuan = chnMat*VrfQuan1;

    % ZF Digital Beamformer
    Vd = chnEqv'*(chnEqv*chnEqv')^-1;
    VdQuan = chnEqvQuan'*(chnEqvQuan*chnEqvQuan')^-1;

    % For Each SNR (Assume the Noise Power equals 1)
    capRef33 = zeros(numUser,snrLen);
    capRef33Quan = zeros(numUser,snrLen);
    for iSnr = 1:snrLen
        % Use Water-Filling Algorithm to Meet Power Constraint
        % (The total beamformer power constraint is ||Vrf*Vd||_F^2=P)
        powAlloc = water_filling_MIMO_ZF(snrVec(iSnr),ones(numUser,1),Vrf*Vd);
        VdNew = Vd*diag(sqrt(powAlloc));

        powAllocQuan = water_filling_MIMO_ZF(snrVec(iSnr),ones(numUser,1),VrfQuan1*VdQuan);
        VdQuanNew = VdQuan*diag(sqrt(powAllocQuan));

        % Calculate Multi-User Capacity, i.e. Sum-Rate of All Users
        for iUser = 1:numUser
            chnUser = chnMat(iUser,:);
            
            powerRx = abs(chnUser*Vrf*VdNew).^2;
            capRef33(iUser,iSnr) = log2( 1 + powerRx(iUser) / ...
                ( 1 + sum(powerRx([1:iUser-1,iUser+1:numUser])) ) );

            powerRxQuan1 = abs(chnUser*VrfQuan1*VdQuanNew).^2;
            capRef33Quan(iUser,iSnr) = log2( 1 + powerRxQuan1(iUser) / ...
                ( 1 + sum(powerRxQuan1([1:iUser-1,iUser+1:numUser])) ) );
        end
    end

    % the Sum Capacity of All Users
    capRef33Sum = capRef33Sum + sum(capRef33)';
    capRef33QuanSum = capRef33QuanSum + sum(capRef33Quan)';

    %% Reference [32]: Hybrid Ananlog Point to Best Path (MUBS) + Digital ZF
    % Analog Beamformer Pointing to Best Path
    Vrf = zeros(numAntTx, numUser);
    for iUser = 1:numUser
        [~, indMaxGain] = max(abs(chnAll(iUser).pathGain));
        Vrf(:,iUser) = 1/sqrt(numAntTx)*exp(1i*2*pi*1/2*(0:numAntTx-1)'*sin(chnAll(iUser).angleDep(indMaxGain)));
    end
    VrfQuan = quantizeByCodebook(codeBook,Vrf);

    % the Equivalent Channel
    chnEqv = chnMat*Vrf;
    chnEqvQuan = chnMat*VrfQuan;

    % ZF Digital Beamformer
    Vd = chnEqv'*(chnEqv*chnEqv')^-1;
    VdQuan = chnEqvQuan'*(chnEqvQuan*chnEqvQuan')^-1;

    % For Each SNR (Assume the Noise Power equals 1)
    capRef32 = zeros(numUser,snrLen);
    capRef32Quan = zeros(numUser,snrLen);
    for iSnr = 1:snrLen
        % Use Water-Filling Algorithm to Meet Power Constraint
        % (The total beamformer power constraint is ||Vrf*Vd||_F^2=P)
        powAlloc = water_filling_MIMO_ZF(snrVec(iSnr),ones(numUser,1),Vrf*Vd);
        VdNew = Vd*diag(sqrt(powAlloc));

        powAllocQuan = water_filling_MIMO_ZF(snrVec(iSnr),ones(numUser,1),VrfQuan*VdQuan);
        VdQuanNew = VdQuan*diag(sqrt(powAllocQuan));

        % Calculate Multi-User Capacity, i.e. Sum-Rate of All Users
        for iUser = 1:numUser
            chnUser = chnMat(iUser,:);

            % Infinite Resolution
            powerRx = abs(chnUser*Vrf*VdNew).^2;
            capRef32(iUser,iSnr) = log2( 1 + powerRx(iUser) / ...
                ( 1 + sum(powerRx([1:iUser-1,iUser+1:numUser])) ) );

            % Finite Resolution
            powerRxQuan = abs(chnUser*VrfQuan*VdQuanNew).^2;
            capRef32Quan(iUser,iSnr) = log2( 1 + powerRxQuan(iUser) / ...
                ( 1 + sum(powerRxQuan([1:iUser-1,iUser+1:numUser])) ) );
        end
    end

    % the Sum Capacity of All Users
    capRef32Sum = capRef32Sum + sum(capRef32)';
    capRef32QuanSum = capRef32QuanSum + sum(capRef32Quan)';

    %% Print the Progress
    disp([num2str(iMC),'-th Monte-Carlo iteration finished!'])

end

%% Averge Data After Monte-Carlo
capFdSum = capFdSum/numMC;
capPropSum = capPropSum/numMC;
capPropQuan1Sum = capPropQuan1Sum/numMC;
capPropQuan2Sum = capPropQuan2Sum/numMC;
capRef33Sum = capRef33Sum/numMC;
capRef33QuanSum = capRef33QuanSum/numMC;
capRef32Sum = capRef32Sum/numMC;
capRef32QuanSum = capRef32QuanSum/numMC;

%% Plot Figure
figure;
plot(snrVecDb,capFdSum,'--k',snrVecDb,capPropSum,'--ok',...
    snrVecDb,capPropQuan1Sum,'--db',snrVecDb,capPropQuan2Sum,'--squareb',...
    snrVecDb,capRef33Sum,'-^r', snrVecDb,capRef33QuanSum,'-vr',...
    snrVecDb,capRef32Sum,'-.xg', snrVecDb,capRef32QuanSum,'-.|g',LineWidth=1);
hold on;
grid on;
xlabel('SNR (dB)');
ylabel('Sum Rate (bits/s/Hz)');
legend("Fully-Digital ZF","Proposed Algorithm for N^{RF}=5, b=\infty",...
    "Proposed Algorithm for N^{RF}=5, b=1, quantize final",...
    "Proposed Algorithm for N^{RF}=5, b=1, quantize intermediate",...
    "Hybrid beamforming in [33], N^{RF}=4, b=\infty", ...
    "Hybrid beamforming in [33], N^{RF}=4, b=1", ...
    "Hybrid Beamforming in [32], N^{RF}=4, b=\infty", ...
    "Hybrid Beamforming in [32], N^{RF}=4, b=1");
