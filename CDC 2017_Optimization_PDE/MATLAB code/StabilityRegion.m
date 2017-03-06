% 2017-02-21
% Stability region as a function of (1) controller parameter: beta, 
%                                   (2) actuator location za
clear;clc;

%% independent variables
% h axis
hstep = 0.01;
h_period = 0:hstep:0.3;
NumH = numel(h_period);

% z axis
BC_ini = 0;                     % boundary condition 1
BC_fin = pi;                    % boundary condition 2
Zstep = pi/100;                 % step of z
z = BC_ini: Zstep: BC_fin;
NumZ = numel(z);      

indep_var = 3;                  % 1: h 
                                % 2: h and beta
                                % 3: h and za
                                
if indep_var == 1
    z_temp = pi/2;
    beta_temp = 2;
    SR = zeros(1, NumH);
elseif indep_var == 2
    z_temp = pi/2;
    beta = 0:0.01:4;
    NumBeta = numel(beta);
    SR = zeros(NumBeta, NumH);
elseif indep_var == 3
    beta_temp = 2;
    SR = zeros(NumZ, NumH);
end



% 
% % z = pi/2;
% NumZ = numel(z);                % number of z intervals




% controller parameter
% beta_temp = 2;
% beta = 0:0.001:2;
% NumBeta = numel(beta);

%% stability function init
% SR = zeros(NumZ, NumH);
% SR = zeros(NumBeta, NumH);
%% Liptchiz constants
% Antsakli's parameters =====
% Kf = 3.6503;
% Kdelta = 0.3155;
% Km = 0.0010;
% Kfm = Kf + Km;
% alpha_c = 1.5;
% beta = 0.9;
% ===========================
alpha_c = 1;
% beta = 2;

%% stability region
Lip_coef = 1;
if indep_var == 1  % F1(h)
    [Kf, Km, Kdelta] = Lipschitz_Constants(z_temp, beta_temp, z, Zstep);
    Kf = Lip_coef*Kf;
    Km = Lip_coef*Km;
    Kdelta = Lip_coef*Kdelta;
    Kfm = Kf + Km;
    
   for hh = 1:NumH
       h_temp = h_period(hh);
       SR(hh) = 1 - alpha_c*(exp(-beta_temp*h_temp) + Kdelta/(Kfm + beta_temp)*(exp(Kfm*h_temp) - exp(-beta_temp*h_temp)));
   end
   
elseif indep_var == 2 % F1(h, beta)
    for bb = 1:NumBeta
        beta_temp = beta(bb);
        [Kf, Km, Kdelta] = Lipschitz_Constants(z_temp, beta_temp, z, Zstep);
        Kf = Lip_coef*Kf;
        Km = Lip_coef*Km;
        Kdelta = Lip_coef*Kdelta;
        Kfm = Kf + Km;

       for hh = 1:NumH
           h_temp = h_period(hh);
           SR(bb,hh) = 1 - alpha_c*(exp(-beta_temp*h_temp) + Kdelta/(Kfm + beta_temp)*(exp(Kfm*h_temp) - exp(-beta_temp*h_temp)));
       end
    end
 
elseif indep_var == 3 % F1(h, za)
    for zz = 1:NumZ
        z_temp = z(zz);
        [Kf, Km, Kdelta] = Lipschitz_Constants(z_temp, beta_temp, z, Zstep);
        Kf = Lip_coef*Kf;
        Km = Lip_coef*Km;
        Kdelta = Lip_coef*Kdelta;
        Kfm = Kf + Km;

       for hh = 1:NumH
           h_temp = h_period(hh);
           SR(zz,hh) = 1 - alpha_c*(exp(-beta_temp*h_temp) + Kdelta/(Kfm + beta_temp)*(exp(Kfm*h_temp) - exp(-beta_temp*h_temp)));
       end
    end
end
   
%%
% temperature profile

if indep_var == 1  % F1(h)
    figure()
    plot(h_period, SR);hold on;
    plot(h_period, zeros(1,NumH));
    xlabel('h','FontSize',30);
    title('$F_1(h)$','Interpreter','latex','FontSize',30);
    set(gca,'FontSize',24);
elseif indep_var == 2  % F1(h, beta)
    figure()
    mesh(h_period,beta,SR);
    xlim([0 max(h_period)]);
    ylim([0 max(beta)]);
    xlabel('h','FontSize',30);
    ylabel('$\beta$','Interpreter','latex','FontSize',30);
    zlabel('$F_1(h,\beta)$','Interpreter','latex','FontSize',30);
    set(gca,'FontSize',24);
    
    figure()
    [c,h] = contourf(h_period,beta,SR); grid on;
    clabel(c,h)
    xlabel('h','FontSize',30);
    ylabel('$\beta$','Interpreter','latex','FontSize',30);
    set(gca,'FontSize',24);
elseif indep_var == 3 % F1(h, za)
    figure()
    mesh(h_period,z,SR);
    xlim([0 max(h_period)]);
    ylim([0 max(z)]);
    xlabel('h','FontSize',30);
    ylabel('$z_a$','Interpreter','latex','FontSize',30);
    zlabel('$F_1(h,z_a)$','Interpreter','latex','FontSize',30);
    set(gca,'FontSize',24);
    
    figure()
    [c,h] = contourf(h_period,z,SR, 'LevelStep', 0.2); grid on;
    clabel(c,h)
    xlabel('h','FontSize',30);
    ylabel('$z_a$','Interpreter','latex','FontSize',30);
    set(gca,'FontSize',24);
end
% mesh(h_period,z,SR);
% % plot(h_period, SR);hold on;
% % plot(h_period, zeros(1,NumH));
% % xlim([0 max(h_period)]);
% % ylim([0 pi]);
% xlabel('h','FontSize',30);
% ylabel('$z_a$','Interpreter','latex','FontSize',30);
% zlabel('$F_1(h,z_a)$','Interpreter','latex','FontSize',30);
% % title('$\max(F_1(h)) = 0.3038$ $a\in [-0.5, 0.5]$','Interpreter','latex','FontSize',30)
% set(gca,'FontSize',24);
% %%
% figure()
% % v = [1,1];
% [c,h] = contourf(h_period,z,SR); grid on;
% clabel(c,h)
% xlabel('h','FontSize',30);
% ylabel('$z_a$','Interpreter','latex','FontSize',30);
% set(gca,'FontSize',24);

