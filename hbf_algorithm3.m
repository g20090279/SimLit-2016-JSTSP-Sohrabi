% last Revised: Apr. 15, 2023.
function Vrf = hbf_algorithm3(chn, numAnt, numRf)

% Initialize the Analog Part of the Hybrid Precoder
Vrf = exp(1i*2*pi*rand(numAnt,numRf));
diff = inf;
eps = 1e-3;
counter = 0;
f_old = 0;

% Iteration
while diff > eps
    for rfInd = 1:numRf
        
        % Use Sherman-Morrison Formula to Separate A rank-k Matrix Inverse
        % as A Sum of A rank-(k-1) Matrix Inverse and An Extra Component
        % Containing rank-1 Matrix.
        Vrfbar = Vrf;
        Vrfbar(:,rfInd) = [];
        Aj = chn*(Vrfbar*Vrfbar')*chn';
        
        for antInd = 1:numAnt
            % This Algorithm Needs numRf>=K+1 in Sherman-Morrison Formula.
            % Otherwise, Aj is not invertible.
            Bj = chn'*Aj^(-2)*chn;
            
            vrf = Vrf(:,rfInd);
            vrf(antInd) = [];

            Bjbar = Bj;
            Bjbar(antInd,:) = [];
            Bjbar(:,antInd) = [];

            ZetaB = Bj(antInd,antInd) + vrf'*Bjbar*vrf; % should be real-valued
            ZetaB = abs(ZetaB); % make sure real-valued by removing computing error

            Bj_row = Bj(antInd,:);
            Bj_row(antInd) = [];
            vrf_col = Vrf(:,rfInd);
            vrf_col(antInd) = [];

            etaB = Bj_row * vrf_col;

            Dj = chn'*Aj^(-1)*chn;
            Djbar = Dj;
            Djbar(antInd,:) = [];
            Djbar(:,antInd) = [];

            ZetaD = Dj(antInd,antInd) + vrf'*Djbar*vrf; % should be real-valued
            ZetaD = abs(ZetaD);

            Dj_row = Dj(antInd,:);
            Dj_row(antInd) = [];

            etaD = Dj_row * vrf_col; 

            cij = (1+ZetaD)*etaB - ZetaB*etaD;
            zij = imag(2*conj(etaB)*etaD);

            if real(cij) >= 0
                phiij = asin(imag(cij)/abs(cij));
            else
                phiij = pi - asin(imag(cij)/abs(cij));
            end

            theta1 = -phiij + asin(zij/abs(cij));
            theta2 = pi - phiij - asin(zij/abs(cij));

            % update
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

    % check convergence
    counter = counter + 1;
    f_new = abs( trace( (chn*(Vrf*Vrf')*chn')^(-1) ));
    if f_old == 0
        f_old = f_new;
        continue;
    else
        diff = f_old - f_new;
        f_old = f_new;
    end
end

% print information
disp(['Converge in ',num2str(counter),' iterations.']);

end