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
ActPos = zeros(1,3);
ActPos(1,1) = pi/3;
ActPos(1,1) = 2*pi/3;
% ActPos(1,3) = 0.3236;
ActPos(1,3) = 0.0377;

SenPos = zeros(1,3);
SenPos(1,1) = pi/4;
SenPos(1,2) = 2*pi/4;
SenPos(1,3) = 3*pi/4;

% ET = 0.2422;                    % model error threshold coefficient
ET = 0.0505;

NumEv = 3;
NumEvs = 3;
NumAc = 3;

% A and Amodel
A = [1 0 0;
     0 -2 0;
     0 0 -5];
 
deltaA = [0.02 0 0;
          0 0.02 0;
          0 0 0.02];
      
Amodel = A + deltaA;
deltaB = 0.5*deltaA;

% C 
Q = zeros(NumEv,NumAc);
for n = 1: 3
    for i = 1:3
    Q(n,i) = sqrt(2/pi)*sin(n*SenPos(1,i)); 
    end
end
C = Q';

% B 
B = zeros(NumEv,NumAc);
for n = 1: NumEv
    for i = 1:NumAc
        B(n,i) = 2*sqrt(2/pi)*sin(n*ActPos(1,i)); 
    end
end
Bmodel = B + deltaB;

% find feedback gain matrix
P = eye(3);
poles = [-100 -100 -100];
K = place(C*Amodel*C^(-1),C*Bmodel,poles);

%% definition of state variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(NumZ, NumT);              % states
Y = zeros(3, NumT);                 % outputs
a = zeros(3, NumT);                 % eigenmodes
U = zeros(3, NumT);                 % inputs

Xmodel = zeros(NumZ, NumT);         % states
Ymodel = zeros(3, NumT);        % outputs
amodel = zeros(3, NumT);       % eigenmodes

% eigenfunctions
Phi = zeros(NumEv,NumZ);
for n = 1: NumEv
    for m = 1:NumZ
       Phi(n,m) = sqrt(2/pi)*sin(n*z(m));
    end
end

% eigenfunctions of model
Phis = zeros(NumEv,NumZ);
for n = 1: NumEvs
    for m = 1: NumZ
        Phis(n,m) = sqrt(2/pi)*sin(n*z(m));
    end
end

ModelError = zeros(NumAc,NumT);     % model estimation error
UF = zeros(NumAc,NumT);             % udpate flag
count = zeros(1,NumAc);             % count update times
NetworkLoad = zeros(1,NumT);        
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PLANT OPERATION PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m: z-axis index
% n: eigenvalue index
% t: time-axis index

%% initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a(1,1) = 0.7;                   % the 1st eigenmode at time=0 is 0.7
amodel(:,1) = a(1:NumEvs,1);

X(:,1) = Phi'*a(:,1);
Y(:,1) = C*a(:,1);
Xmodel(:,1) = Phis'*amodel(:,1);
Ymodel(:,1) = C*amodel(:,1);
U(:,1) = K*Ymodel(:,1);
%% the big loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 2: NumT
    dYmodel = C*Amodel*C^(-1)*Ymodel(:,t-1) - C*B*U(:,t-1);
    Ymodel(:,t) = Ymodel(:,t-1) + dYmodel*Tstep;
    
    dY = C*A*C^(-1)*Y(:,t-1) - C*B*U(:,t-1);
    Y(:,t) = Y(:,t-1) + dY*Tstep;
    a(:,t) = C^(-1)*Y(:,t);
    X(:,t) = Phis'*a(:,t);
    
    for i = 1:3
       ModelError(i,t) = norm(Y(i,t) - Ymodel(i,t)); 
       if ModelError(i,t) > ET*norm(Y(i,t))    % update yhat i
           Ymodel(i,t) = Y(i,t);
           UF(i,t) = 1;
           count(1,i) = count(1,i) + 1;
       end
    end
    NetworkLoad(1,t) = UF(1,t) + UF(2,t) + UF(3,t);
    
    amodel(:,t) = C^(-1)*Ymodel(:,t);
    Xmodel(:,t) = Phis'*amodel(:,t);
    U(:,t) = K*Ymodel(:,t);

end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PLOT SIMULATION RESULTES %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% temperature profile
figure()
% mesh(T,z,Xmodel);
mesh(T,z,X);
axis([0 tFinal 0 pi -0.6 0.6])
xlabel('Time,t','FontSize',30);
ylabel('z','FontSize',30)
zlabel('','Interpreter','latex','String','$\bar x(t,z)$','FontSize',30)
set(gca,'FontSize',24);

% control actions
figure()
plot(T,U(1,:));hold on;
plot(T,U(2,:));hold on;
plot(T,U(3,:));
legend('U1','U2','U3','FontSize',24)
xlabel('Time,t','FontSize',30);
ylabel('U','FontSize',30)
set(gca,'FontSize',24);

% % updates
% figure()
% plot(T,UF(1,:)); hold on;
% plot(T,UF(2,:)); hold on;
% plot(T,UF(3,:))
% legend('Update1','Update2','Update3','FontSize',24)
% xlabel('Time, t','FontSize',30);
% ylabel('Update','FontSize',30)
% set(gca,'FontSize',24);
% 
% % network load
% figure()
% plot(T,NetworkLoad); 
% axis([0 tFinal 0 3.2]);
% xlabel('Time, t','FontSize',30);
% ylabel('Network load','FontSize',30)
% set(gca,'FontSize',24);
% 
% % model estimation errors
% figure()
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