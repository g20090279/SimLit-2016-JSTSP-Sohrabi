% Test basic linear algebra
% Last revised: Mar 22, 2023
% Zekai Liang, liang@ice.rwth-aachen.de
%%
Nt = 16;
Nr = 4;
Nrf = 6;
H = randn(Nr,Nt)+1i*rand(Nr,Nt);
Vrf = 1/sqrt(Nt)*exp(1i*2*pi*rand([Nt,Nrf]));


%%
j = randi(Nrf);  % 1:Nr

Heff = H*Vrf;
A = Heff*Heff';

Vrfj = Vrf;
Vrfj(:,j) = [];

Heffj = H*Vrfj;
Aj = Heffj*Heffj';
eig(Aj)

H1 = H*Vrf(:,j);
Ar = H1*H1';

A2 = Aj + Ar;
% check if A2==A

%%
Ainv = inv(A);  % Nrf>=Nr, otherwise, singular
Ainv*A    % check if A^-1*A = I

Ajinv = inv(Aj); % Nrf>Nr, otherwise, singular
Ajinv*Aj  % check if Aj^-1*Aj = I
A2inv = Ajinv - Ajinv*Ar*Ajinv/(1+Vrf(:,j)'*H'*Ajinv*H*Vrf(:,j));