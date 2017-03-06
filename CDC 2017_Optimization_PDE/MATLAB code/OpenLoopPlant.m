% Simulate open-loop system

clear;clc;

% initialize simulation constants
Tstep = 0.001;          % specifiy step of time
tFinal = 5;             % final time
NumEv = 500;             % number of eigenvalues used
BC_ini = 0;             % boundary condition 1
BC_fin = pi;            % boundary condition 2
Zstep = pi/100;         % step of z
% colocated actuators and sensors
ActPos = zeros(1,3);
ActPos(1,1) = 1;        % position of actuator 1
ActPos(1,2) = 2;        % position of actuator 2
ActPos(1,3) = 3;        % position of actuator 3
SenPos = zeros(1,3);
SenPos(1,1) = 1;        % position of sensor 1
SenPos(1,2) = 2;        % position of sensor 2
SenPos(1,3) = 3;        % position of sensor 3
% parameters
betaT = 200;            % dimensionless heat of rxn
betaU = 2;              % dimensionless heat transfer coefficient
gamma = 4;              % dimensionless activation energy
theta1 = 0;             % parametric uncertainty in the heat of rxn
theta2 = 0;             % the maximum value of theta1
theta2_Pos = 0.125*pi;  % the location of the expected uncertainty
alpha = 1;              % difussion coefficient

% initial simulation parameters
% construct time axis
T = 0: Tstep: tFinal;
NumT = numel(T);        % number of time intervals
% construct z axis
z = BC_ini: Zstep: BC_fin;
NumZ = numel(z);        % number of z intervals
X = zeros(NumZ, NumT);  % state
Y = zeros(3, NumT);
a = zeros(NumEv, NumT); % eigenmodes
% aUS = zeros(NumEvs, NumT);  % unstable eigenmodes

% eigenvalues: lambda = -betaU - n^2
lambda = zeros(1, NumEv);
for n = 1: NumEv
   lambda(1,n) = -betaU - n^2;
end

% state-space expression
% m: z-axis index
% n: eigenvalue index
Ainf = diag(lambda);
Binf = zeros(NumEv,3);
Qinf = zeros(NumEv,3);
Winf = zeros(NumEv,1);

% Binf matrix
for n = 1: NumEv
   for i = 1:3
      Binf(n,i) = 2*sqrt(2/pi)*sin(n*ActPos(1,i)); 
   end
end

% Qinf matrix
for n = 1: NumEv
   for i = 1:3
      Qinf(n,i) = sqrt(2/pi)*sin(n*SenPos(1,i)); 
   end
end

% Winf matrix
for n = 1: NumEv
   Winf(n,1) = betaU*sqrt(2/pi)*sin(n*theta2_Pos);
end

% initial condition 
% m: z-axis index
% n: eigenvalue index
a(1,1) = 0.7;           % the 1st eigenmode at time=0 is 0.7
for m = 1: NumZ
   for n = 1: NumEv
      X(m,1) = X(m,1) + a(n,1)*sqrt(2/pi)*sin(n*z(m)); 
   end
end
for i = 1:3
    for n = 1: NumEv
        Y(i,n) = Y(i,n) + a(n,1)*sqrt(2/pi)*sin(n*SenPos(1,i)); 
    end
end
% eigenfunctions
Phi = zeros(NumEv,NumZ);
for n = 1: NumEv
    for m = 1:NumZ
       Phi(n,m) = sqrt(2/pi)*sin(n*z(m));
    end
end
% big loop
% t: time-axis index
for t = 2: NumT
    % eigenmodel a
    dan = NonlinearFunction(Ainf,a(:,t-1),gamma,betaT,Phi,X(:,t-1),Zstep);
    a(:,t) = a(:,t-1) + Tstep*dan;
    X(:,t)=Phi'*a(:,t);         % The actual state follows x(z)=sum a(z)*phi(z)
    Y(:,t)=Qinf'*a(:,t);
end

% temperature profile
mesh(T,z,X);
xlabel('Time, t','FontSize',30);
ylabel('z','FontSize',30)
zlabel('Temperature','FontSize',30);
set(gca,'FontSize',24);







