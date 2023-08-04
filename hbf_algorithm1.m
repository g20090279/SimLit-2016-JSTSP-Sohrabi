function [V,C] = hbf_algorithm1(F,alpha,numColV)
% HBF_ALGORITHM1 solves the following optimization problem by an iterative
% coordinate descent algorithm.
%
%                   max_Vrf    C = log2(det( I + alpha*V'*F*V ))
%                    s.t.      |V(i,j)|^2 <= 1, for all i, j
%
% where V is a numRowV-by-numColV matrix, F is a numColV-by-numColV matrix.
% alpha is a scalar containing SNR value in linear.
% The objective function can be rewritten to a simpler form by decomposing
% matrix F blockwise, i.e.
%
%                          | scalar,  | row vector |
%                     F =  | ---------|------------|
%                          |  column  |            |
%                          |  vector, |    matrix  |
%
% Input(s):
% - F: a square matrix.
% - alpha: a scalar.
% - numRowV: the number of rows for matrix V
%
% Output(s):
% - V: a matrix.
% - C: a scalar.

% Make Sure F is A Square Matrix
numRowV = size(F);
if numRowV(1)~=numRowV(2), error('F should be a square matrix!'); end
numRowV = numRowV(1);

% Convergence Threshold
epsilon = 10^-3;

% Initialize V
V = exp(1i*2*pi*rand([numRowV,numColV]));

% Calculate Initial Capacity
C = abs(log2(det( eye(numColV) + alpha*V'*F*V )));

% Initialize Variables
diff = inf;
count = 0;

while diff > epsilon
    for iCol = 1:numColV
        VrfBar = V;
        VrfBar(:,iCol) = [];
        Cj = eye(numColV-1) + alpha*VrfBar'*F*VrfBar;
        Gj = alpha*F - alpha^2*F*VrfBar*Cj^(-1)*VrfBar'*F;
        for iRow = 1:numRowV
            % Obtain Vector Gj(i,l) and Vrf(l,j) Where l~=i
            gj = Gj(iRow,:);
            gj(iRow) = [];
            vrfj = V(:,iCol);
            vrfj(iRow) = [];

            % Eq.(14) Determines Vrf(i,j) value
            V(iRow,iCol) = exp(1i*angle(gj*vrfj));
        end
    end

    % Calculate the New Capacity, Assuming Vd*Vd'=I
    cNew = abs(log2(det( eye(numColV) + alpha*V'*F*V )));
    diff = abs(cNew-C)/C;
    C = cNew;

    % Count How Many Iterations to Converge
    count = count+1;
end

% Display Which Iteration We Are In
% disp(['For alpha=',num2str(alpha),': converged after ',num2str(count),' iterations.']);
end