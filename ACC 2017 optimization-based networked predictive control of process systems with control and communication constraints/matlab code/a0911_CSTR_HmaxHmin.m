%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2016-09-07 for 2017 ACC
% Characterize Hmax and Hmin for 2D CSTR system
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

deltaA = A - Ahat;
deltaB = B - Bhat;

k1 = 4.2e4:10:5e4;
k2 = 1400:5:3200;

NumK1 = length(k1);
NumK2 = length(k2);

indkk = ones(NumK1,NumK2);
H_temp = zeros(NumK1,NumK2,4);
Hmax = zeros(NumK1,NumK2);
Hmin = zeros(NumK1,NumK2);

H = 0:0.005:0.3;
NumH = length(H);

Wx = eye(2);
Wu = eye(1);
K_lqr = lqr(Ahat, Bhat, Wx, Wu); % 1.0e+04 *[4.4118    0.2000]
%%

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
           lambda(p,q,h) = max(abs(eig(M)));
            
%            if abs(lambda(p,q,h)-1) < 0.2
%                H_temp(p,q,indkk(p,q)) = H(h);
%                indkk(p,q) = indkk(p,q) + 1;
%            end
       end
    end
    waitbar(p/NumK1,hwaitbar,strcat(num2str(p/NumK1*100),'%'));
end
%%
for p = 1:NumK1
   for q = 1:NumK2
       for h = 1:NumH-1
          if lambda(p,q,h)>=1 && lambda(p,q,h+1)<=1
             Hmin_temp = H(h);
             if Hmin(p,q) < Hmin_temp
                Hmin(p,q) = Hmin_temp;
             end
          end
          if lambda(p,q,h)<=1 && lambda(p,q,h+1)>1
             Hmax_temp = H(h);
             if Hmax(p,q) < Hmax_temp
                Hmax(p,q) = Hmax_temp;
             end
          end
       end
   end
end

for p = 1:NumK1
   for q = 1:NumK2
      if Hmax(p,q) == 0
         Hmax(p,q) = NaN; 
      end
      if Hmin(p,q) == 0
         Hmin(p,q) = NaN; 
      end
   end
end
%%
figure()
[c,ch] = contourf(k1,k2,Hmax'); grid on;
clabel(c,ch)
xlabel('$k_1$','Interpreter','LaTex'); 
ylabel('$k_2$','Interpreter','LaTex'); 
zlabel('$H_{\max}$', 'Interpreter','Latex')

figure()
mesh(k1,k2,Hmax');
xlabel('$k_1$','Interpreter','LaTex'); 
ylabel('$k_2$','Interpreter','LaTex'); 
zlabel('$H_{\max}$', 'Interpreter','Latex')

figure()
[c,ch] = contourf(k1,k2,Hmin'); grid on;
clabel(c,ch)
xlabel('$k_1$','Interpreter','LaTex'); 
ylabel('$k_2$','Interpreter','LaTex'); 
zlabel('$H_{\min}$', 'Interpreter','Latex')

figure()
mesh(k1,k2,Hmin');
xlabel('$k_1$','Interpreter','LaTex'); 
ylabel('$k_2$','Interpreter','LaTex'); 
zlabel('$H_{\min}$', 'Interpreter','Latex')

figure()
mesh(k1,k2,Hmax');hold on;
mesh(k1,k2,Hmin');
xlabel('$k_1$','Interpreter','LaTex'); 
ylabel('$k_2$','Interpreter','LaTex'); 
zlabel('$h$', 'Interpreter','Latex')

toc



