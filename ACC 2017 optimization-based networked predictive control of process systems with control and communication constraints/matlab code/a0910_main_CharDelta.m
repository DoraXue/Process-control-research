%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2016-09-10 for 2017 ACC
% Characterize Delta contour from K for 2D system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;
tic;

% % tunable parameters

% lqr: scan k1
Wx = eye(2);
Wu = eye(1);
% K_lqr = lqr(Ahat, Bhat, Wx, Wu); % 2.3762    0.0000
k1 = 4.2e4:5;
k2 = -5:0.01:10;
H = 0:0.01:3;
NumK1 = length(k1);
NumK2 = length(k2);

Hmax = zeros(NumK1, NumK2);
Hmin = zeros(NumK1, NumK2);
%%
fontsize=20;
set(0,'defaultaxesfontsize',fontsize); % 0 is for current MATLAB session
set(0,'defaulttextfontsize',fontsize);

% plant parameters
A = diag([1,-1]);
B = [1; 1];

% model parameters
Ahat = diag([0.8,-0.8]);
Bhat = [1.1; 1.2];

deltaA = A - Ahat;
deltaB = B - Bhat;

NumH = length(H);

lambda = zeros(NumK1, NumK2, NumH);
Is = [ones(size(A)) zeros(size(A)); zeros(size(A)) zeros(size(A))];

hwaitbar = waitbar(0,strcat(num2str(0/NumK1*100),'%'));

for p = 1:NumK1
    for q = 1:NumK2
       for h = 1:NumH
           K = [k1(p) k2(q)];
           Lambda = [(A - B*K) B*K;
                    (deltaA - deltaB*K) (Ahat + deltaB*K)];
           M = Is*expm(Lambda*H(h))*Is;
           lambda(p,q,h) = max(norm(eig(M)));
           
           if abs(lambda(p,q,h)-1) < 0.1
               if Hmax(p,q) == 0
                  Hmax(p,q) = H(h);
               else
                  Hmax(p,q) = max(Hmax(p,q),H(h)); 
               end
               
               if Hmin(p,q) == 0
                  Hmin(p,q) = H(h);
               else
                  Hmin(p,q) = min(Hmin(p,q),H(h)); 
               end
           end
       end
    end
   waitbar(p/NumK1,hwaitbar,strcat(num2str(p/NumK1*100),'%'));
end
%%
figure()
[c,h] = contourf(k1,k2,Hmax'); grid on;
clabel(c,h)
xlabel('$k_1$','Interpreter','LaTex'); 
ylabel('$k_2$','Interpreter','LaTex'); 
zlabel('$H_{max}$', 'Interpreter','Latex')

figure()
mesh(k1,k2,Hmax');
axis([min(k1) max(k1) min(k2) max(k2) 0 3])
xlabel('$k_1$','Interpreter','LaTex'); 
ylabel('$k_2$','Interpreter','LaTex'); 
zlabel('$H_{min}$', 'Interpreter','Latex')


figure()
[c,h] = contourf(k1,k2,Hmin'); grid on;
clabel(c,h)
xlabel('$k_1$','Interpreter','LaTex'); 
ylabel('$k_2$','Interpreter','LaTex'); 
zlabel('$H_{min}$', 'Interpreter','Latex')

figure()
mesh(k1,k2,Hmin');
axis([min(k1) max(k1) min(k2) max(k2) 0 3])
xlabel('$k_1$','Interpreter','LaTex'); 
ylabel('$k_2$','Interpreter','LaTex'); 
zlabel('$H_{min}$', 'Interpreter','Latex')






