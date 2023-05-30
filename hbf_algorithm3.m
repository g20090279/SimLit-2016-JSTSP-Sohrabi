% Algorithm 3 (for multiuser)
% last Revised: May 25, 2023.
% Zekai Liang, liang@ice.rwth-aachen.de
% Iterate between power P and analog beamformer Vrf

function [Vrf,Vd,P] = hbf_algorithm3(chn, numAnt, numRf, numUser, snr)

% Initialize the Analog Part of the Hybrid Precoder
P = eye(numUser);  % initial power is an identity matrix
chnTilde = P^(-1/2)*chn;
Vrf = exp(1i*2*pi*rand(numAnt,numRf));

% Initialize Variables for Vrf and P Iterative Optimization
counter = 0;
diff = inf;
eps = 1e-3;
fOld = 0;

% Optimize Original Power Constraint Without Relaxation
while diff > eps

    % Optimize Vrf with A Fixed Power P
    while diff > eps

        for rfInd = 1:numRf

            % Use Sherman-Morrison Formula to Separate A rank-k Matrix Inverse
            % as A Sum of A rank-(k-1) Matrix Inverse and An Extra Component
            % Containing rank-1 Matrix.
            Vrfbar = Vrf;  % Vrfbar is the submatrix of Vrf with j(rfInd)-th column removed
            Vrfbar(:,rfInd) = [];
            Aj = chnTilde*(Vrfbar*Vrfbar')*chnTilde';

            for antInd = 1:numAnt
                % This Algorithm Needs numRf>=K+1 in Sherman-Morrison Formula.
                % Otherwise, Aj is not invertible.
                Bj = chnTilde'*Aj^(-2)*chnTilde;

                vrfBar = Vrf(:,rfInd);
                vrfBar(antInd) = [];

                Bjbar = Bj;
                Bjbar(antInd,:) = [];
                Bjbar(:,antInd) = [];

                ZetaB = Bj(antInd,antInd) + vrfBar'*Bjbar*vrfBar; % should be real-valued
                ZetaB = abs(ZetaB); % make sure real-valued by removing computing error

                BjRow = Bj(antInd,:);
                BjRow(antInd) = [];
                vrfCol = Vrf(:,rfInd);
                vrfCol(antInd) = [];

                etaB = BjRow * vrfCol;

                Dj = chnTilde'*Aj^(-1)*chnTilde;
                DjBar = Dj;
                DjBar(antInd,:) = [];
                DjBar(:,antInd) = [];

                ZetaD = Dj(antInd,antInd) + vrfBar'*DjBar*vrfBar; % should be real-valued
                ZetaD = abs(ZetaD);

                DjRow = Dj(antInd,:);
                DjRow(antInd) = [];

                etaD = DjRow * vrfCol;

                cij = (1+ZetaD)*etaB - ZetaB*etaD;
                zij = imag(2*conj(etaB)*etaD);

                if real(cij) >= 0
                    phiij = asin(imag(cij)/abs(cij));
                else
                    phiij = pi - asin(imag(cij)/abs(cij));
                end

                theta1 = -phiij + asin(zij/abs(cij));
                theta2 = pi - phiij - asin(zij/abs(cij));

                % Update
                opt1 = exp(-1i*theta1);
                opt2 = exp(-1i*theta2);
                f1 = ZetaB+2*real(conj(opt1)*etaB)/(1+ZetaD+2*real(conj(opt1)*etaD));
                f2 = ZetaB+2*real(conj(opt2)*etaB)/(1+ZetaD+2*real(conj(opt2)*etaD));
                if f1>=f2
                    Vrf(antInd,rfInd) = opt1;
                else
                    Vrf(antInd,rfInd) = opt2;
                end
            end
        end

        % Check Convergence of Analog Precoder Vrf
        counter = counter + 1;
        fNew = abs( numAnt*trace( (chnTilde*(Vrf*Vrf')*chnTilde')^-1 ) );  % Careful! P is not invertible!
        if fOld == 0
            fOld = fNew;
            continue;
        else
            diff = fOld - fNew;
            fOld = fNew;
        end
    end

    % After Convergence of Vrf, Caculate Digital Precoder with Zero-Forcing
    Vd = Vrf'*chn'*(chn*(Vrf*Vrf')*chn')^(-1);

    % Water-Filling Algorithm for Power Allocation
    P = diag( water_filling_MIMO_ZF(snr,ones(numUser,1),Vrf*Vd) );
    
    % Note That P Has To Be Invertible After Water-Filling. Otherwise, The
    % Following Algorithm Doesn't Work.
    if ~all(diag(P)>0)
        break;
    end

    chnTilde = P^(-1/2)*chn;
    fNew = abs( numAnt*trace( (chnTilde*(Vrf*Vrf')*chnTilde')^-1 ));  % careful !
    if fOld == 0
        fOld = fNew;
        continue;
    else
        diff = fOld - fNew;
        fOld = fNew;
    end

end

% End of Function, Print Information
P = diag(P);
disp(['Converge in ',num2str(counter),' iterations.']);

end