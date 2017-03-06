clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% information of coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct time axis
Tstep = 1e-3;                   % specifiy step of time
tFinal = 5;                     % final time
T = 0: Tstep: tFinal;
NumT = numel(T);                % number of time intervals

% construct z axis
BC_ini = 0;                     % boundary condition 1
BC_fin = pi;                    % boundary condition 2
Zstep = pi/100;                 % step of z
z = BC_ini: Zstep: BC_fin;
NumZ = numel(z);                % number of z intervals

%% colocated actuators and sensors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumAc = 3;
ActPos = zeros(1,NumAc);
for ac = 1: NumAc
   ActPos(1,ac) = ac;
%    ActPos(1,ac) = ac*pi/(NumAc+1); 
end
SenPos = ActPos;                % colocated sensors and actuators

%% information of the plant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaTreal = 200;                % dimensionless heat of rxn
betaU = 2;                      % dimensionless heat transfer coefficient
gamma = 4;                      % dimensionless activation energy
theta1 = 0;                     % parametric uncertainty in the heat of rxn
theta2 = 1;                     % the maximum value of theta1
theta2_Pos = [0.125*pi 0.625*pi 0.825*pi];          % the location of the expected uncertainty
NumDis = numel(theta2_Pos);     % number of point disturbances
thetab = theta2;
alpha = 1;                      % difussion coefficient
betaT = betaTreal + theta1;     % betaT with uncertainty used in models

%% information of eigenvalues %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumEv = NumAc;                  % number of eigenvalues used for infinite dimensional plant
NumEvs = NumAc;                 % number of unstable eigenvalues/order of reduced-order system

% eigenvalues: lambda = -betaU - n^2
lambda = zeros(1, NumEv);
for n = 1: NumEv
   lambda(1,n) = -betaU - alpha*n^2;
end

%% definition of state variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% under Broadcast-based communication strategy
X_B = zeros(NumZ, NumT);          % states
Y_B = zeros(NumAc, NumT);         % outputs
a_B = zeros(NumEv, NumT);         % eigenmodes
U_B = zeros(NumAc, NumT);         % inputs

% under Adaptive communication strategy
X_A = zeros(NumZ, NumT);          % states
Y_A = zeros(NumAc, NumT);         % outputs
a_A = zeros(NumEv, NumT);         % eigenmodes
U_A = zeros(NumAc, NumT);         % inputs

% state-space expression
% m: z-axis index
% n: eigenvalue index
Ainf = diag(lambda);
Binf = zeros(NumEv,NumAc);
Qinf = zeros(NumEv,NumAc);
Winf = zeros(NumEv,NumDis);

% Binf matrix
for n = 1: NumEv
   for i = 1:NumAc
      Binf(n,i) = 2*sqrt(2/pi)*sin(n*ActPos(1,i)); 
   end
end

% Qinf matrix
for n = 1: NumEv
   for i = 1:NumAc
      Qinf(n,i) = sqrt(2/pi)*sin(n*SenPos(1,i)); 
   end
end

% Winf matrix
for n = 1: NumEv
    for d = 1: NumDis
        Winf(n,d) = betaU*sqrt(2/pi)*sin(n*theta2_Pos(d));
    end
end

% eigenfunctions
Phi = zeros(NumEv,NumZ);
for n = 1: NumEv
    for m = 1:NumZ
       Phi(n,m) = sqrt(2/pi)*sin(n*z(m));
    end
end

%% initialization of models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% under Broadcast-based communication strategy
Xmodel_B = zeros(NumZ, NumT);         % states
Ymodel_B = zeros(NumAc, NumT);        % outputs
amodel_B = zeros(NumEvs, NumT);       % eigenmodes
Umodel_B = zeros(NumAc, NumT);        % inputs

% under Adaptive communication strategy
Xmodel_A = zeros(NumZ, NumT);         % states
Ymodel_A = zeros(NumAc, NumT);        % outputs
amodel_A = zeros(NumEvs, NumT);       % eigenmodes
Umodel_A = zeros(NumAc, NumT);        % inputs

Amodel = zeros(NumEvs,NumEvs);
Bmodel = zeros(NumEvs,NumAc);
Qmodel = zeros(NumEvs,NumAc);
Wmodel = zeros(NumEvs,1);           % a zeros matrix for models
Phis = zeros(NumEvs,NumZ);

deltaA = [0.02 0 0;
          0 0.02 0;
          0 0 0.02];
Amodel = diag(lambda(1,1:NumEvs)) + deltaA;

% Bmodel matrix
for n = 1:NumEvs
   for i = 1:NumAc
      Bmodel(n,i) = 2*sqrt(2/pi)*sin(n*ActPos(1,i)); 
   end
end

% Qmodel matrix
for n = 1: NumEvs
   for i = 1:NumAc
      Qmodel(n,i) = sqrt(2/pi)*sin(n*SenPos(1,i)); 
   end
end

% eigenfunctions
for n = 1: NumEvs
    for m = 1: NumZ
        Phis(n,m) = sqrt(2/pi)*sin(n*z(m));
    end
end

%% controller parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1*eye(NumEvs);              % quadratic Lyapunov function
L = 1000;                        % observer gain
ModelError = zeros(NumAc,NumT); % model estimation error
ET = [0.8 0.8 0.8];             % model error thresholds
dV_A = zeros(1,NumT);
UF_B = zeros(NumAc,NumT);       % udpate flag using Broadcast-based communication
UF_A = zeros(1,NumT);       % udpate flag using Adaptive communication
count_B = zeros(1,NumAc);         % count update times, 1st row: B; 2nd row: A
count_A = 0;
Net = zeros(1,NumT);
TotalData_B = 0;
TotalData_A = 0;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PLANT OPERATION PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m: z-axis index
% n: eigenvalue index
% t: time-axis index

%% initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_B(1,1) = 0.7;                   % the 1st eigenmode at time=0 is 0.7
a_A(1,1) = 0.7;
amodel_B(:,1) = a_B(1:NumEvs,1);
amodel_A(:,1) = a_A(1:NumEvs,1);

X_B(:,1) = Phi'*a_B(:,1);
Y_B(:,1) = Qinf'*a_B(:,1);
Xmodel_B(:,1) = Phis'*amodel_B(:,1);
Ymodel_B(:,1) = Qmodel'*amodel_B(:,1);

X_A(:,1) = Phi'*a_A(:,1);
Y_A(:,1) = Qinf'*a_A(:,1);
Xmodel_A(:,1) = Phis'*amodel_A(:,1);
Ymodel_A(:,1) = Qmodel'*amodel_A(:,1);

norm_B = zeros(1,NumT);
norm_A = zeros(1,NumT);
V_A(1,:) = transpose(a_A(:,1))*P*a_A(:,1);

%% the big loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 2: NumT
    % under Broad-cast based communication strategy %%%%%%%%%%%%%%%%%%%%%%%
    % find control action
    NonFunObserver_B = NonlinearFunction(Amodel,amodel_B(:,t-1),gamma,betaT,Phis,Xmodel_B(:,t-1),Zstep);
    dVdx1_B = 2*transpose(amodel_B(:,t-1))*P;
    LfV_B = dVdx1_B*(NonFunObserver_B + L*(Ymodel_B(:,t-1) - Qmodel*amodel_B(:,t-1)));
    LgV_B = dVdx1_B*Bmodel;
    % dV = LfV + LgV*U(:,t-1); 
    LwV_B = dVdx1_B*Phis*(exp(-gamma./(1.+Xmodel_B(:,t-1))) - exp(-gamma))*Zstep;
    U_B(:,t) = NonlinearController(amodel_B(:,t-1),LfV_B,LgV_B,LwV_B,thetab);
    
    % state and model state evolution
    da_B = NonlinearFunction(Ainf,a_B(:,t-1),gamma,betaTreal,Phi,X_B(:,t-1),Zstep);
    a_B(:,t) = a_B(:,t-1) + Tstep*(da_B + Binf*U_B(:,t));
    X_B(:,t) = Phi'*a_B(:,t);
    Y_B(:,t) = Qinf'*a_B(:,t);
    
    damodel_B = NonlinearFunction(Amodel,amodel_B(:,t-1),gamma,betaT,Phis,Xmodel_B(:,t-1),Zstep);
    amodel_B(:,t) = amodel_B(:,t-1) + Tstep*(damodel_B + Bmodel*U_B(:,t));
    Xmodel_B(:,t) = Phis'*amodel_B(:,t);
    Ymodel_B(:,t) = Qmodel'*amodel_B(:,t);
    
    % update models if necessary
    % calculate model estimation errors
    for n = 1: NumAc
       ModelError(n,t) = norm(Ymodel_B(n,t) - Y_B(n,t));
    end
    % update y1
    if ModelError(1,t) > ET(1)*norm(Y_B(1,t)) % update Ymodel(:,t)
        count_B(1,1) = count_B(1,1) + 1;
        UF_B(1,t) = 1;
        Ymodel_B(1,t) = Y_B(1,t);
    end
    % update y2
    if ModelError(2,t) > ET(2)*norm(Y_B(2,t)) % update Ymodel(:,t)
        count_B(1,2) = count_B(1,2) + 1;
        UF_B(2,t) = 1;
        Ymodel_B(2,t) = Y_B(2,t);
    end
    % update y3
    if ModelError(3,t) > ET(3)*norm(Y_B(3,t)) % update Ymodel(:,t)
        count_B(1,3) = count_B(1,3) + 1;
        UF_B(3,t) = 1;
        Ymodel_B(3,t) = Y_B(3,t);
    end
        amodel_B(:,t) = Qmodel^(-1)*Ymodel_B(:,t);
        Xmodel_B(:,t) = Phis'*amodel_B(:,t);  
        Net(1,t) = UF_B(1,t) + UF_B(2,t) + UF_B(2,t);
   norm_B(1,t) = norm(Y_B(:,t));
   if t*Tstep >= 2
      TotalData_B = TotalData_B + Net(1,t);
   end
    
    
    % under Adaptive communication strategy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find control action
    NonFunObserver_A = NonlinearFunction(Amodel,amodel_A(:,t-1),gamma,betaT,Phis,Xmodel_A(:,t-1),Zstep);
    dVdx1_A = 2*transpose(amodel_A(:,t-1))*P;
    LfV_A = dVdx1_A*(NonFunObserver_A + L*(Ymodel_A(:,t-1) - Qmodel*amodel_A(:,t-1)));
    LgV_A = dVdx1_A*Bmodel;
    dV_A(1,t) = LfV_A + LgV_A*U_A(:,t-1); 
    LwV_A = dVdx1_A*Phis*(exp(-gamma./(1.+Xmodel_A(:,t-1))) - exp(-gamma))*Zstep;
    U_A(:,t) = NonlinearController(amodel_A(:,t-1),LfV_A,LgV_A,LwV_A,thetab);
    
    % state evolution
    da_A = NonlinearFunction(Ainf,a_A(:,t-1),gamma,betaTreal,Phi,X_A(:,t-1),Zstep);
    a_A(:,t) = a_A(:,t-1) + Tstep*(da_A + Binf*U_A(:,t));
    X_A(:,t) = Phi'*a_A(:,t);
    Y_A(:,t) = Qinf'*a_A(:,t);
%     V_A(:,t) = transpose(a_A(:,t))*P*a_A(:,t);
%     dV_A = V_A(:,t) - V_A(:,t-1);
    
    
    % update models if necessary
    % calculate the overall Vdot
    norm_A(1,t) = norm(Y_A(:,t));
    if norm_A(1,t) >= 0.01
        if dV_A(1,t) >= 0       % if outside the terminal set and dotV >= 0
            count_A = count_A + 1;
            UF_A(1,t) = 1;
            Ymodel_A(:,t) = Y_A(:,t);
            amodel_A(:,t) = Qmodel^(-1)*Ymodel_A(:,t); 
            Xmodel_A(:,t) = Phis'*amodel_A(:,t); 
        end
    else                             % do not update any y
        damodel_A = NonlinearFunction(Amodel,amodel_A(:,t-1),gamma,betaT,Phis,Xmodel_A(:,t-1),Zstep);
        amodel_A(:,t) = amodel_A(:,t-1) + Tstep*(damodel_A + Bmodel*U_A(:,t));
        Xmodel_A(:,t) = Phis'*amodel_A(:,t);
        Ymodel_A(:,t) = Qmodel'*amodel_A(:,t);
    end
    if t*Tstep >= 2
      TotalData_A = TotalData_A + 3*UF_A(1,t);
   end
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PLOT SIMULATION RESULTES %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% temperature profile
figure()
% mesh(T,z,Xmodel);
mesh(T,z,X_B);
axis([0 5 0 pi -0.2 0.6])
xlabel('Time,t','FontSize',30);
ylabel('z','FontSize',30)
zlabel('','Interpreter','latex','String','$\bar x(t,z)$','FontSize',30)
set(gca,'FontSize',24);

% control actions
figure()
plot(T,U_B(1,:));hold on;
plot(T,U_B(2,:));hold on;
plot(T,U_B(3,:));
legend('U1','U2','U3','FontSize',24)
xlabel('Time,t','FontSize',30);
ylabel('U','FontSize',30)
set(gca,'FontSize',24);

% updates using Broad-cast based communication strategy
figure()
plot(T,UF_B(1,:)); hold on;
plot(T,UF_B(2,:)); hold on;
plot(T,UF_B(3,:))
set(gca,'YLim',[0 1]);%X????????
set(gca,'YTick',[0:1:1]);%?????????
% set(gca,'XTickLabel',[0:0.1:1.5]);%?????? 
legend('Update1','Update2','Update3','FontSize',24)
xlabel('Time,t','FontSize',30);
ylabel('Update','FontSize',30)
set(gca,'FontSize',24);

% updates using Adaptive communication strategy
figure()
plot(T,UF_A)
set(gca,'YLim',[0 1]);%X????????
set(gca,'YTick',[0:1:1]);%?????????
% set(gca,'XTickLabel',[0:0.1:1.5]);%?????? 
% legend('Update1','Update2','Update3','FontSize',24)
xlabel('Time,t','FontSize',30);
ylabel('Update','FontSize',30)
set(gca,'FontSize',24);

figure()
% plot(T,3*UF_A)
plot(T,Net)
set(gca,'YLim',[-0.5 3.5]);%X????????
set(gca,'YTick',[0:1:3]);%?????????
% set(gca,'XTickLabel',[0:0.1:1.5]);%?????? 
% legend('Update1','Update2','Update3','FontSize',24)
xlabel('Time,t','FontSize',30);
ylabel('Network Load','FontSize',30)
set(gca,'FontSize',24);

figure()
plot(T,norm_B)
xlabel('Time,t','FontSize',30);
ylabel('||y||','FontSize',30)
set(gca,'FontSize',24);

% % model estimation errors
% figure(4)
% subplot(3,1,1);
% plot(T,ModelError(1,:)); 
% % xlabel('Time, t','FontSize',30);
% % ylabel('Model estimation errors','FontSize',30)
% set(gca,'FontSize',24);
% 
% subplot(3,1,2);
% plot(T,ModelError(2,:)); 
% % xlabel('Time, t','FontSize',30);
% ylabel('Model estimation errors','FontSize',30)
% set(gca,'FontSize',24);
% 
% subplot(3,1,3);
% plot(T,ModelError(3,:)); 
% % legend('Error1','Error2','Error3','FontSize',24)
% xlabel('Time, t','FontSize',30);
% % ylabel('Model estimation errors','FontSize',30)
% set(gca,'FontSize',24);