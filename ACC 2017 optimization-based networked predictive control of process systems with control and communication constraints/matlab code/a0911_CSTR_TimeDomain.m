%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2016-09-09 for 2017 ACC
% Time-domain: optimization-based controller test
% 2D CSTR
% used for paper draft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;
tic;

fontsize=20;
set(0,'defaultaxesfontsize',fontsize); % 0 is for current MATLAB session
set(0,'defaulttextfontsize',fontsize);

% tunable parameters
horizon = 1;
Wx = eye(2);
Wu = 1;
Wh = 1e3;

coef = 0.7;

dist = 3; 
dis_int = 4;
dis_end = dis_int + 0.01;

% plant parameters
A = [ -1.7398462744754825032212573211112, -0.027234749658550006368548841124702;
      147.96925489509650064425146422224,    4.4469499317100012737097682249403];
B = [0; 0.0042];

Ahat = [-1.8 0;
        150 5];
Bhat = [0; 0.005];

deltaA = A - Ahat;
deltaB = B - Bhat;
%%
load('0911_k1.mat');plot
load('0911_k2.mat');
% load('0907_lambda.mat');
load('0911_Hmax.mat');
load('0911_Hmin.mat');

NumK1 = length(k1);
NumK2 = length(k2);

% simulation parameters
Tstep = 0.01;
time = 0:Tstep:8;
NumT = length(time);

x = zeros(2, NumT);
x(1,1) = 0.1;
xm = x;
K_act = zeros(2, NumT);

K_LQR = lqr(Ahat, Bhat, Wx, Wu);
K_lqr = [K_LQR(1)*ones(1,NumT); K_LQR(2)*ones(1,NumT)];

%%
update = zeros(1,NumT);
update(1) = 1;

u = zeros(1,NumT);

CostJ = zeros(1, NumT);
CostJ(1) = x(:,1)'*Wx*x(:,1) + Wu*u(1)^2 + Wh/1;

Delta = 0;

hwaitbar = waitbar(0,strcat(num2str(0/NumT*100),'%'));
t_flag = 2;

period = zeros(1,NumT);
% period(1) = 0;
%%
for t = 2:NumT
    if t == floor(t_flag + (Delta/Tstep)) % solve optimization
        % update model state
        xm(:,t-1) = x(:,t-1);
        update(t) = 1;
        t_flag = t;
        
        J = zeros(NumK1,NumK2);
        Jmin_temp = 1e20;
        
        for indk1 = 1:NumK1
            for indk2 = 1:NumK2
                if ~isnan(Hmin(indk1,indk2)) && ~isnan(Hmax(indk1,indk2))
                   K_temp = [k1(indk1) k2(indk2)];
                   Delta_temp = Hmin(indk1,indk2) + coef*(Hmax(indk1,indk2) - Hmin(indk1,indk2));

                   J(indk1,indk2) = 0;
                   xm0 = xm(:,t-1);
                   xm_pre = xm0;
                   for indt = 1:(horizon/Tstep)
                       dxmdt = (Ahat - Bhat*K_temp)*xm_pre;
                       xm_now = xm_pre + dxmdt*Tstep;

                       u_temp = K_temp*xm_pre;

                       J(indk1,indk2) = J(indk1,indk2) + (xm_now'*Wx*xm_now + Wu*u_temp^2)*Tstep;

                       xm_pre = xm_now;
                   end
                   J(indk1,indk2) = J(indk1,indk2) + Wh/Delta_temp;
                   
                   if J(indk1,indk2) < Jmin_temp
                       Jmin_temp = J(indk1,indk2);
                       k_ind1 = indk1;
                       k_ind2 = indk2;
                   end
                end
%                 waitbar((t*indk1*indk2)/(NumT*NumK1*NumK2),hwaitbar,strcat(num2str((t*indk1*indk2)/(NumT*NumK1*NumK2)*100),'%'));

            end
        end
        

        K_act(:,t) = [k1(k_ind1) k2(k_ind2)];
        
        Delta = Hmin(k_ind1,k_ind2) + coef*(Hmax(k_ind1,k_ind2) - Hmin(k_ind1,k_ind2));
        period(1,t) = Delta/Tstep;
    else
        K_act(:,t) = K_act(:,t-1);
        period(1,t) = period(1,t-1);
    end
    
    % introduce disturbance
    if t >= dis_int/Tstep && t <= dis_end/Tstep
        dxdt = A*x(:,t-1) - B*K_act(:,t)'*xm(:,t-1) + dist;
    else
        dxdt = A*x(:,t-1) - B*K_act(:,t)'*xm(:,t-1);
    end
    u(1,t) = -K_act(:,t)'*xm(:,t-1);
%     dxdt = A*x(:,t-1) - B*K_act(:,t)'*xm(:,t-1);
    x(:,t) = x(:,t-1) + dxdt*Tstep;
    
    dxmdt = Ahat*xm(:,t-1) - Bhat*K_act(:,t)'*xm(:,t-1);
    xm(:,t) = xm(:,t-1) + dxmdt*Tstep;
    
    CostJ(t) = x(:,t)'*Wx*x(:,t) + Wu*u(t)^2 + Wh/Delta;

    
    waitbar(t/NumT,hwaitbar,strcat(num2str(t/NumT*100),'%'));
end

UF_opt = sum(update)/NumT*100
%%
x_lqr = zeros(2,NumT);
x_lqr(1) = 0.1;
xm_lqr = x_lqr;
update_lqr = zeros(1,NumT);
update_lqr(1) = 1;
u_lqr = zeros(1,NumT);

Delta_lqr = 0.17 + coef*(0.21 - 0.17);

CostJ_lqr = zeros(1, NumT);
CostJ_lqr(1) = x_lqr(:,1)'*Wx*x_lqr(:,1) + Wu*u_lqr(1)^2 + Wh/1;


for t = 2:NumT
    if mod(t,floor(Delta_lqr/Tstep)) == 0 % solve optimization
        xm_lqr(:,t-1) = x_lqr(:,t-1);   
        update_lqr(:,t) = 1;
    end
    
    % introduce disturbance
    if t >= dis_int/Tstep && t <= dis_end/Tstep
        dxdt = A*x_lqr(:,t-1) - B*K_lqr(:,t)'*xm_lqr(:,t-1) + dist;
    else
        dxdt = A*x_lqr(:,t-1) - B*K_lqr(:,t)'*xm_lqr(:,t-1);
    end
    u_lqr(1,t) =  - K_lqr(:,t)'*xm_lqr(:,t-1);
%     dxdt = A*x_lqr(:,t-1) - B*K_lqr(:,t)'*xm_lqr(:,t-1);
    x_lqr(:,t) = x_lqr(:,t-1) + dxdt*Tstep;
    
    
    dxmdt = Ahat*xm_lqr(:,t-1) - Bhat*K_lqr(:,t)'*xm_lqr(:,t-1);
    xm_lqr(:,t) = xm_lqr(:,t-1) + dxmdt*Tstep;
    CostJ_lqr(t) = x_lqr(:,t)'*Wx*x_lqr(:,t) + Wu*u_lqr(t)^2 + Wh/Delta_lqr;

end
UF_lqr = sum(update_lqr)/NumT*100

CA = x(1,:) + 0.575;
T = x(2,:) + 395.05;
%%
figure()
subplot(2,1,1)
plot(time, x(1,:),'LineWidth',2); hold on;
plot(time, x_lqr(1,:),'r--');
legend('optimal','LQR') 
ylabel('x');

subplot(2,1,2)
plot(time, x(2,:),'LineWidth',2); hold on;
plot(time, x_lqr(2,:),'r--');
legend('optimal','LQR')
xlabel('time'); ylabel('x');
%%
figure()
subplot(2,1,1)
plot(time, CA,'LineWidth',2); hold on;
plot(time, x_lqr(1,:)+0.575,'r--');
legend('optimal','LQR') 
ylabel('$C_A$(mol/L)','Interpreter','LaTex');

subplot(2,1,2)
plot(time, T,'LineWidth',2); hold on;
plot(time, x_lqr(2,:)+395.05,'r--');
legend('optimal','LQR')
xlabel('$t$(min)','Interpreter','LaTex'); 
ylabel('$T$(K)','Interpreter','LaTex');
%%
figure()
subplot(2,1,1);
stem(time, update,'LineWidth',2,'Marker','none'); 
axis([0 max(time) -0.1 1.1])
ylabel('optimal');

subplot(2,1,2);
stem(time, update_lqr,'LineWidth',2,'Marker','none');
axis([0 max(time) -0.1 1.1])
xlabel('$t$(min)','Interpreter','LaTex'); 
ylabel('LQR');
%%
figure()
subplot(2,1,1)
plot(time, K_act(1,:),'LineWidth',2);hold on;
plot(time, K_lqr(1,:), 'r--' ,'LineWidth',2)
legend('optimal', 'LQR')
ylabel('k_1');

subplot(2,1,2)
plot(time, K_act(2,:),'LineWidth',2);hold on;
plot(time, K_lqr(2,:), 'r--' ,'LineWidth',2)
legend('optimal', 'LQR')
xlabel('$t$(min)','Interpreter','LaTex'); 
ylabel('k_2');

%%
figure()
plot(time,u,'LineWidth',2);hold on;
plot(time, u_lqr, 'r')
legend('optimal', 'LQR')
xlabel('$t$(min)','Interpreter','LaTex'); 
ylabel('u');

toc
%%
figure()
plot(time,period*Tstep,'LineWidth',2);hold on;
plot(time, Delta_lqr*ones(1,NumT), 'r' ,'LineWidth',2)
legend('optimal', 'LQR')
xlabel('$t$(min)','Interpreter','LaTex'); 
ylabel('$h$(min)','Interpreter','LaTex');
