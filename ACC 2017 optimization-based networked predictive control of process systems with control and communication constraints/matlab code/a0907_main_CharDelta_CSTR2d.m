%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2016-09-07 for 2017 ACC
% Characterize Delta from K for 2D system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize=20;
set(0,'defaultaxesfontsize',fontsize); % 0 is for current MATLAB session
set(0,'defaulttextfontsize',fontsize);

% model on nael's book p.90
clear;clc;close all;
tic;
A = [ -1.7398462744754825032212573211112, -0.027234749658550006368548841124702;
      147.96925489509650064425146422224,    4.4469499317100012737097682249403];
B = [0; 0.0042];

Ahat = [-1.8 0;
        150 5];
Bhat = [0; 0.005];

% k1 = 4e4:10:1.2e5;
% k2 = 0.2e4;

% k1 = 4.4118e4;
k1 = 4.4e4;       % for boundary check
% k1 = 4.41e4;       % for boundary check
k2 = 1400:3200;

% k2 = 1400;          % for boundary check
% k2 = 2000;          % for boundary check
% k1 = 4.2e4:10:5e4;

H = 0:0.01:0.3;

Wx = eye(2);
Wu = eye(1);
K_lqr = lqr(Ahat, Bhat, Wx, Wu); % 1.0e+04 *[4.4118    0.2000]
%%
NumK = length(k1);

deltaA = A - Ahat;
deltaB = B - Bhat;

NumH = length(H);

lambda = zeros(NumK, NumH);
Is = [ones(size(A)) zeros(size(A)); zeros(size(A)) zeros(size(A))];

hwaitbar = waitbar(0,strcat(num2str(0/NumK*100),'%'));

for p = 1:NumK
   for h = 1:NumH
       K = [k1(p) k2];
%        K = [k1 k2(p)];
       Lambda = [(A - B*K) B*K;
                (deltaA - deltaB*K) (Ahat + deltaB*K)];
       M = Is*expm(Lambda*H(h))*Is;
       lambda(p,h) = max(abs(eig(M)));

   end
   waitbar(p/NumK,hwaitbar,strcat(num2str(p/NumK*100),'%'));
end
%%
figure()
v = [1,1];
[c,h] = contourf(k1,H,lambda',v); grid on;
% axis([4.1e4 7e4 0 0.26])
clabel(c,h)
xlabel('$k_1$','Interpreter','LaTex'); 
ylabel('$h$(min)','Interpreter','LaTex'); 
% title('H for LQR');
% zlabel('$\lambda$', 'Interpreter','Latex')
%%
% figure()
% mesh(k1,H,lambda');
% % axis([min(k2) max(k2) min(H) max(H) 0 110])
% xlabel('$k_1$','Interpreter','LaTex'); 
% ylabel('$h$(min)','Interpreter','LaTex'); 
% zlabel('$\lambda_{max}$', 'Interpreter','Latex')
% 
% 
% 
% 



