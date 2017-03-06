% find a constant K that stabilizes the closed-loop system regardless of Za
% Decision: K = 25
clear;clc;

% construct z axis
BC_ini = 0;                     % boundary condition 1
BC_fin = pi;                    % boundary condition 2
Zstep = pi/100;                 % step of z
z = BC_ini: Zstep: BC_fin;

mu_b = 0;                     % half-width of finite support

NumZ = numel(z);                % number of z intervals

% colocated actuators and sensors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumAc = 1;

% information of the plant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaTreal = 80;                 % dimensionless heat of rxn
betaU = 2;                      % dimensionless heat transfer coefficient
betaUmodel = 2.05;
gamma = 4;                      % dimensionless activation energy
theta1 = 5;                     % parametric uncertainty in the heat of rxn
alpha = 1;                      % difussion coefficient
betaT = betaTreal + theta1;     % betaT with uncertainty used in models

% information of eigenvalues %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumEv = NumAc;                  % number of eigenvalues used for infinite dimensional plant

% eigenvalues: lambda = -betaU - n^2
lambda = -betaU - alpha;

Ainf = lambda;
Amodel = -betaUmodel - alpha;
% Binf = betaU*sqrt(2/pi)*sin(ActPos);

% eigenfunctions
Phi = zeros(NumEv,NumZ);
for n = 1: NumEv
    for m = 1:NumZ
       Phi(n,m) = sqrt(2/pi)*sin(n*z(m));
    end
end

Phi_model = zeros(NumEv,NumZ);
for n = 1: NumEv
    for m = 1:NumZ
       Phi_model(n,m) = sqrt(2/pi)*sin(n*z(m));
    end
end

BC_ini = 0;                     % boundary condition 1
BC_fin = pi;                    % boundary condition 2
Zstep = pi/100;                 % step of z
z = BC_ini: Zstep: BC_fin;
NumZ = length(z);

pp = -5;
K = 25;
CL_pole = zeros(1,NumZ);

for i = 2:NumZ-1
    ActPos = z(i);
    X = Phi'*ActPos;
    da(i) = NonlinearFunction(Ainf,ActPos,gamma,betaTreal,Phi,X,Zstep);
    
    damodel = NonlinearFunction(Amodel,ActPos,gamma,betaT,Phi_model,X,Zstep);
    deltaf(i) = da(i) - damodel;
    
    U = (-da(i) - beta*a(i))/Bmodel;
    ga(i) =  Binf*U;
    gamodel(i) = -da(i) - beta*a(i);
    deltag(i) = ga(i) - gamodel(i);
    
    
    
    Binf = betaU*mu_b*sqrt(1/(2*pi))*(cos(ActPos - mu_b) - cos(ActPos + mu_b));
    Bmodel = betaUmodel*sqrt(2/pi)*sin(ActPos);

    K(i) = place(Amodel, Bmodel,pp);
    CL_pole(i) = Ainf - Binf*K(i);
end

figure()
plot(z, K); hold on;
xlim([0 pi]);
xlabel('z_a','FontSize',30);
ylabel('K','FontSize',30)
title('Find constant K that works for all z_a');
set(gca,'FontSize',24);

figure()
plot(z, CL_pole); hold on;
xlim([0 pi]);
xlabel('z_a','FontSize',30);
ylabel('CL pole','FontSize',30)
title('Closed-loop pole location');
set(gca,'FontSize',24);


