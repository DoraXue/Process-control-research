clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% information of coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct time axis
Tstep = 1e-3;           % specifiy step of time
tFinal = 5;             % final time
T = 0: Tstep: tFinal;
NumT = numel(T);        % number of time intervals

% construct z axis
BC_ini = 0;             % boundary condition 1
BC_fin = pi;            % boundary condition 2
Zstep = pi/100;         % step of z
z = BC_ini: Zstep: BC_fin;
NumZ = numel(z);        % number of z intervals

% colocated actuators and sensors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ActPos = zeros(1,3);
% ActPos(1,1) = 1;        % position of actuator 1
% ActPos(1,2) = 2;        % position of actuator 2
% ActPos(1,3) = 3;        % position of actuator 3
NumAc = 4;
ActPos = zeros(1,NumAc);
for ac = 1: NumAc
%     ActPos(1,ac) = ac;
   ActPos(1,ac) = ac*pi/(NumAc+1); 
end
SenPos = ActPos;        % colocated sensors and actuators

% information of the plant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaTreal = 100;            % dimensionless heat of rxn
betaU = 2;              % dimensionless heat transfer coefficient
gamma = 4;              % dimensionless activation energy
theta1 = 2;             % parametric uncertainty in the heat of rxn
theta2 = 2;             % the maximum value of theta1
theta2_Pos = 0.125*pi;  % the location of the expected uncertainty
thetab = 2;
alpha = 1;              % difussion coefficient
betaT = betaTreal + theta1; % betaT with uncertainty used in models

% information of eigenvalues %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumEv = 50;             % number of eigenvalues used for infinite dimensional plant
NumEvs = NumAc;             % number of unstable eigenvalues/order of reduced-order system

% eigenvalues: lambda = -betaU - n^2
lambda = zeros(1, NumEv);
for n = 1: NumEv
   lambda(1,n) = -betaU - n^2;
end

% definition of state variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(NumZ, NumT);  % states
Y = zeros(NumAc, NumT);     % outputs
a = zeros(NumEv, NumT); % eigenmodes
U = zeros(NumAc, NumT);     % inputs

% state-space expression
% m: z-axis index
% n: eigenvalue index
Ainf = diag(lambda);
Binf = zeros(NumEv,NumAc);
Qinf = zeros(NumEv,NumAc);
Winf = zeros(NumEv,1);

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
   Winf(n,1) = betaU*sqrt(2/pi)*sin(n*theta2_Pos);
end

% eigenfunctions
Phi = zeros(NumEv,NumZ);
for n = 1: NumEv
    for m = 1:NumZ
       Phi(n,m) = sqrt(2/pi)*sin(n*z(m));
    end
end

% initialization of models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xmodel = zeros(NumZ, NumT);     % states
Ymodel = zeros(NumAc, NumT);        % outputs
amodel = zeros(NumEvs, NumT);   % eigenmodes
Umodel = zeros(NumAc, NumT);        % inputs

Amodel = zeros(NumEvs,NumEvs);
Bmodel = zeros(NumEvs,NumAc);
Qmodel = zeros(NumEvs,NumAc);
Wmodel = zeros(NumEvs,1);
Phis = zeros(NumEvs,NumZ);

Amodel = diag(lambda(1,1:NumEvs));

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

% Wmodel matrix
for n = 1: NumEv
   Wmodel(n,1) = betaU*sqrt(2/pi)*sin(n*theta2_Pos);
end

% eigenfunctions
for n = 1: NumEvs
    for m = 1: NumZ
        Phis(n,m) = sqrt(2/pi)*sin(n*z(m));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PLANT OPERATION PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m: z-axis index
% n: eigenvalue index
% t: time-axis index

P = 1*eye(NumEvs);             % quadratic Lyapunov function
L = 100;                % observer gain

% initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a(1,1) = 0.7;           % the 1st eigenmode at time=0 is 0.7
amodel(:,1) = a(1:NumEvs,1);

X(:,1) = Phi'*a(:,1);
Y(:,1) = Qinf'*a(:,1);
Xmodel(:,1) = Phis'*amodel(:,1);
Ymodel(:,1) = Qmodel'*amodel(:,1);

% the big loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 2: NumT
    % find control action
    NonFunObserver = NonlinearFunction(Amodel,amodel(:,t-1),gamma,betaT,Phis,Xmodel(:,t-1),Zstep);
    dVdx1 = 2*transpose(amodel(:,t-1))*P;
    LfV = dVdx1*(NonFunObserver + L*(Ymodel(:,t-1) - Qmodel*amodel(:,t-1)));
    LgV = dVdx1*Bmodel;
    % dV = LfV + LgV*U(:,t-1); 
    LwV = dVdx1*Phis*(exp(-gamma./(1.+Xmodel(:,t-1))) - exp(-gamma))*Zstep;
    U(:,t) = NonlinearController(amodel(:,t-1),LfV,LgV,LwV,thetab);
    
    % state and model state evolution
    da = NonlinearFunction(Ainf,a(:,t-1),gamma,betaTreal,Phi,X(:,t-1),Zstep);
    a(:,t) = a(:,t-1) + Tstep*(da + Binf*U(:,t));
    X(:,t) = Phi'*a(:,t);
    Y(:,t) = Qinf'*a(:,t);
    
    damodel = NonlinearFunction(Amodel,amodel(:,t-1),gamma,betaT,Phis,Xmodel(:,t-1),Zstep);
    amodel(:,t) = amodel(:,t-1) + Tstep*(damodel + Bmodel*U(:,t));
    Xmodel(:,t) = Phis'*amodel(:,t);
    Ymodel(:,t) = Qmodel'*amodel(:,t);
    
    % update models
    Ymodel(:,t) = Y(:,t);
    amodel(:,t) = Qmodel^(-1)*Y(:,t);
    Xmodel(:,t) = Phis'*amodel(:,t);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PLOT SIMULATION RESULTES %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% temperature profile
figure(1)
% mesh(T,z,Xmodel);
mesh(T,z,X);
xlabel('Time, t','FontSize',30);
ylabel('z','FontSize',30)
zlabel('Temperature','FontSize',30);
set(gca,'FontSize',24);

% control actions
figure(2)
plot(T,U(1,:));hold on;
plot(T,U(2,:));hold on;
plot(T,U(3,:));
legend('U1','U2','U3','FontSize',24)
xlabel('Time, t','FontSize',30);
ylabel('U','FontSize',30)
set(gca,'FontSize',24);
