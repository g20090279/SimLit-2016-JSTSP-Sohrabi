function powAlloc = water_filling_MIMO_ZF(powMax,powNoise,V)
% WATER_FILLING_MIMO_ZF returns the optimized transmitter power allocation 
% scheme by using water-filling algorithm in a MIMO transmission model. The
% signal model of this MIMO model is
%                        y = H * V * P * x + n,
% where x is a Ns-by-1 transmit signal vector, n and y are noise vector and
% receive signal vector with size  Nr-by-1, respectively. P is a diagonal
% matrix, whose diagonal elements are the optimally allocated transmit
% power for each stream. V is the precoder with size Nt-by-Ns. H is the
% channel matrix with size Nr-by-Nt.
%
% The optimization problem is then formulated
%         max_{p_1,p_2,...,p_Ns}   sum(1+p_k/powNoise_k)
%                       s.t.       Tr(Q*P) <= powMax,
% where Q=V'*V.
%
% The solution is
%              p_k = 1/q_kk * max( 1/lambda - q_kk * powNoise_k , 0 ),
% where q_kk is the k-th diagonal element of Q, and lambda is the water
% line.
%                     __
%               __   |  |         __
%              |  |__|  |      __|  |
%              |  |  |  |__   |  |  |
%            __|  |  |  |  |  |  |  |
%           |  |  |  |  |  |__|  |  |      <----initial water line 1/lambda
%           |  |  |  |  |  |  |  |  |
%     -------------------------------------
%                            ^--- The worst channel with power min(q_kk*powNoise_k)
%
% Inputs:
% - powMax (scalar): The maximum power allowed to be allocated in each 
%     stream.
% - powNoise (scalar): The noise power in linear unit.
% - V (matrix): The beamformer with size Nt-by-Ns. Its power will be 
%      considered during the power allocation. Assume 
%      P = diag(p1,p2,...,pNs) is the allocated power matrix.
%                       Tr(V'*V*P) <= powMax
%      Default: identity matrix.
%
% Outputs:
% - p (col. vector): The allocated power

numStream = size(V,2);

% Predefined Error for Ending Algorithm
eps = 1e-5;

% Acquire Diagonal Elements of the Power
Q = V'*V;  % the diagonal elements must be real-valued
qkk = abs(diag(Q));  % take out the diagonal elements
powEff = powNoise .* qkk;

% Initialize the Waterline to the Lowest One
waterLine = min(powEff);
powAlloc = max( waterLine - powEff, 0 );
powAllocTot = sum(powAlloc);
counter = 0;

% Algorithm: The Idea Is to Raise the Water Level Iteratively to Approach
% the Power Limit.
while abs(powMax - powAllocTot) > eps

    waterLine = waterLine + (powMax-powAllocTot)/numStream;
    powAlloc =  max( waterLine - powEff, 0 );
    powAllocTot = sum(powAlloc);  % power constraint Tr(QP) <= powMax
    counter = counter + 1;

end

disp(['Water-filling finishes in ',num2str(counter),' iterations']);

% The Final Values of Power Allocation
powAlloc = 1./qkk .* powAlloc;

end