clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% information of coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct time axis
Tstep = 1e-3;                   % specifiy step of time
tFinal = 10;                     % final time
T = 0: Tstep: tFinal;
NumT = numel(T);                % number of time intervals

% construct z axis
BC_ini = 0;                     % boundary condition 1
BC_fin = pi;                    % boundary condition 2
Zstep = pi/50;                 % step of z
z = BC_ini: Zstep: BC_fin;
NumZ = numel(z);                % number of z intervals

%% colocated actuators and sensors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumAc = 2;
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
theta1 = 2;                     % parametric uncertainty in the heat of rxn
theta2 = max(theta1);           % the maximum value of theta1
theta2_Pos = 0.125*pi;
% theta2_Pos = [0.125*pi 0.625*pi 0.825*pi];          % the location of the expected uncertainty
NumDis = numel(theta2_Pos);     % number of point disturbances
thetab = theta2;
alpha = 1;                      % difussion coefficient
betaT = betaTreal + theta1;     % betaT with uncertainty used in models

%% information of eigenvalues %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumEv = 2;                  % number of eigenvalues used for infinite dimensional plant
NumEvs = 2;                 % number of unstable eigenvalues/order of reduced-order system

% eigenvalues: lambda = -betaU - n^2
lambda = zeros(1, NumEv);
for n = 1: NumEv
   lambda(1,n) = -betaU - alpha*n^2;
end

%% definition of state variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% under Broadcast-based communication strategy
X = zeros(NumZ, NumT);          % states
Y = zeros(NumAc, NumT);         % outputs
a = zeros(NumEv, NumT);         % eigenmodes
U = zeros(NumAc, NumT);         % inputs

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
Xmodel = zeros(NumZ, NumT);         % states
Ymodel = zeros(NumAc, NumT);        % outputs
amodel = zeros(NumEvs, NumT);       % eigenmodes
Umodel = zeros(NumAc, NumT);        % inputs

Amodel = zeros(NumEvs,NumEvs);
Bmodel = zeros(NumEvs,NumAc);
Qmodel = zeros(NumEvs,NumAc);
Wmodel = zeros(NumEvs,1);           % a zeros matrix for models
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

% eigenfunctions
for n = 1: NumEvs
    for m = 1: NumZ
        Phis(n,m) = sqrt(2/pi)*sin(n*z(m));
    end
end

%% controller parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1*eye(NumEvs);              % quadratic Lyapunov function
L = 1000;                       % observer gain
ModelError = zeros(NumAc,NumT); % model estimation error
ET = [0.8 0.8 0.8];             % model error thresholds
dV = zeros(1,NumT);
UF = zeros(NumAc,NumT);         % udpate flag for each sensor
NetworkLoad = zeros(1,NumT);    % total number of broadcasting sensors
mode = ones(1,NumT);            % 1 = broadcast-based; 2 = adaptive
Terminal = 0.0007;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PLANT OPERATION PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m: z-axis index
% n: eigenvalue index
% t: time-axis index

%% initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a(1,1) = 0.7;                   % the 1st eigenmode at time=0 is 0.7
amodel(:,1) = a(1:NumEvs,1);

X(:,1) = Phi'*a(:,1);
Y(:,1) = Qinf'*a(:,1);
Xmodel(:,1) = Phis'*amodel(:,1);
Ymodel(:,1) = Qmodel'*amodel(:,1);

normA = zeros(1,NumT);

%% the big loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 2: NumT
    % under Broad-cast based communication strategy %%%%%%%%%%%%%%%%%%%%%%%
    % find control action
    NonFunObserver = NonlinearFunction(Amodel,amodel(:,t-1),gamma,betaT,Phis,Xmodel(:,t-1),Zstep);
    dVdx1 = 2*transpose(amodel(:,t-1))*P;
    LfV = dVdx1*(NonFunObserver + L*(Ymodel(:,t-1) - Qmodel*amodel(:,t-1)));
    LgV = dVdx1*Bmodel;
    dV(1,t) = LfV + LgV*U(:,t-1); 
    LwV = dVdx1*Phis*(exp(-gamma./(1.+Xmodel(:,t-1))) - exp(-gamma))*Zstep;
    U(:,t) = NonlinearController(amodel(:,t-1),LfV,LgV,LwV,thetab);
    
    % state and model state evolution
    da = NonlinearFunction(Ainf,a(:,t-1),gamma,betaTreal,Phi,X(:,t-1),Zstep);
    a(:,t) = a(:,t-1) + Tstep*(da + Binf*U(:,t));
    X(:,t) = Phi'*a(:,t);
    Y(:,t) = Qinf'*a(:,t);
    
    % check if any output has entered the terminal set
    for n = 1:NumAc
       if norm(Y(n,t)) <= Terminal
          mode(1,t) = 2; 
       end
    end
    
    if mode(1,t) == 1       % broadcast-based mode
        damodel = NonlinearFunction(Amodel,amodel(:,t-1),gamma,betaT,Phis,Xmodel(:,t-1),Zstep);
        amodel(:,t) = amodel(:,t-1) + Tstep*(damodel + Bmodel*U(:,t));
        Xmodel(:,t) = Phis'*amodel(:,t);
        Ymodel(:,t) = Qmodel'*amodel(:,t);
        
        % update models if necessary
        % calculate model estimation errors
        for n = 1: NumAc
            ModelError(n,t) = norm(Ymodel(n,t) - Y(n,t));
        end
        % update y1
        if ModelError(1,t) > ET(1)*norm(Y(1,t)) % update Ymodel(:,t)
            NetworkLoad(1,t) = NetworkLoad(1,t) + 1;
            UF(1,t) = 1;
            Ymodel(1,t) = Y(1,t);
        end
        % update y2
        if ModelError(2,t) > ET(2)*norm(Y(2,t)) % update Ymodel(:,t)
            NetworkLoad(1,t) = NetworkLoad(1,t) + 1;
            UF(2,t) = 1;
            Ymodel(2,t) = Y(2,t);
        end
%         % update y3
%         if ModelError(3,t) > ET(3)*norm(Y(3,t)) % update Ymodel(:,t)
%             NetworkLoad(1,t) = NetworkLoad(1,t) + 1;
%             UF(3,t) = 1;
%             Ymodel(3,t) = Y(3,t);
%         end
            amodel(:,t) = Qmodel^(-1)*Ymodel(:,t);
            Xmodel(:,t) = Phis'*amodel(:,t);
            
    else        % adaptive mode
        normA(1,t) = norm(Y(:,t));
        if normA(1,t) >= Terminal
            if dV(1,t) >= 0       % if outside the terminal set and dotV >= 0
%                 UF(:,t) = [1;1;1];
                NetworkLoad(1,t) = 3;
                Ymodel(:,t) = Y(:,t);
                amodel(:,t) = Qmodel^(-1)*Ymodel(:,t); 
                Xmodel(:,t) = Phis'*amodel(:,t); 
            end
        else                  % do not update any y
            damodel = NonlinearFunction(Amodel,amodel(:,t-1),gamma,betaT,Phis,Xmodel(:,t-1),Zstep);
            amodel(:,t) = amodel(:,t-1) + Tstep*(damodel + Bmodel*U(:,t));
            Xmodel(:,t) = Phis'*amodel(:,t);
            Ymodel(:,t) = Qmodel'*amodel(:,t);
        end     
%         end
    end
%     Trac(a(1,t),a(2,t),t) = a(3,t);
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PLOT SIMULATION RESULTES %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state trajectory
figure()

r = Terminal;
theta = 0:pi/50:2*pi;
x = r*cos(theta);
y = r*sin(theta);
plot(x,y,'-',0,0'.'); hold on;


plot(Y(1,:),Y(2,:));
axis([-2*Terminal 2*Terminal -2*Terminal 2*Terminal])
xlabel('y1','FontSize',30);
ylabel('y2','FontSize',30)
% zlabel('Temperature','FontSize',30);
set(gca,'FontSize',24);

