
% 2017-02-21
% Stability region as a function of (1) controller parameter: beta, 
%                                   (2) actuator location za
clear;clc;

%% independent variables
% h axis
hstep = 0.001;
h_period = 0:hstep:0.5;
NumH = numel(h_period);

% z axis
BC_ini = 0;                     % boundary condition 1
BC_fin = pi;                    % boundary condition 2
Zstep = pi/100;                 % step of z
z = BC_ini: Zstep: BC_fin;
NumZ = numel(z);      

% beta_temp = 2;
SR = zeros(NumZ, NumH);

alpha_c = 1;
% information of the plant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaTreal = 80;                 % dimensionless heat of rxn
betaU = 2;                      % dimensionless heat transfer coefficient
betaUmodel = 2.05;
gamma = 4;                      % dimensionless activation energy
theta1 = 5;                     % parametric uncertainty in the heat of rxn
alpha = 1;                      % difussion coefficient
betaT = betaTreal + theta1;     % betaT with uncertainty used in models


%% stability region
Lip_coef = 1;
list_Kf = zeros(1,NumZ);
list_Km = zeros(1,NumZ);
list_Kdelta = zeros(1,NumZ);
list_Kfm = zeros(1,NumZ);

K = 15;                         % constant feedback gain
for zz = 1:NumZ
    z_temp = z(zz);
    [Kf, Km, Kdelta] = Lipschitz_Constants_ConstantK(z_temp);
    list_Kf(zz) = Lip_coef*Kf;
    list_Km(zz) = Lip_coef*Km;
    list_Kdelta(zz) = Lip_coef*Kdelta;
    list_Kfm(zz) = Kf + Km;
    
    Amodel = -betaUmodel - alpha;
    Bmodel = betaUmodel*sqrt(2/pi)*sin(z_temp);
    beta_temp = abs(Amodel - Bmodel*K);
    
    for hh = 1:NumH
        h_temp = h_period(hh);
        SR(zz,hh) = 1 - alpha_c*(exp(-beta_temp*h_temp) + ...
            list_Kdelta(zz)/(list_Kfm(zz) + beta_temp)*(exp(list_Kfm(zz)*h_temp) - exp(-beta_temp*h_temp)));
    end
    zz/NumZ
end
%%

figure()
mesh(h_period,z,SR);
xlim([0 max(h_period)]);
ylim([0 max(z)]);
zlim([-1.5 1]);
title(strcat('K= ', num2str(K)));
xlabel('$h$','Interpreter','latex','FontSize',30);
ylabel('$z_a$','Interpreter','latex','FontSize',30);
zlabel('$F_1(h,z_a)$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',24);

figure()
[c,h] = contourf(h_period,z,SR, 'LevelStep', 0.2); grid on;
colormap(hot)
clabel(c,h)
title(strcat('K= ', num2str(K)));
xlabel('$h$','Interpreter','latex','FontSize',30);
ylabel('$z_a$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',24);
