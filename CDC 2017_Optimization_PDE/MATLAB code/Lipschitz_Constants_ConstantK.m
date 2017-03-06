function [Kf, Km, Kdelta] = Lipschitz_Constants_ConstantK(Za)
% Za = pi/2;
% beta = 2;
K = 15;                        % constant feedback gain

% construct z axis
BC_ini = 0;                     % boundary condition 1
BC_fin = pi;                    % boundary condition 2
Zstep = pi/100;                 % step of z
z = BC_ini: Zstep: BC_fin;
% =========================================================================
% 2017-03-02 Introduced uncertainty in b(z) 
% b(z) = mu*sqrt(1/(2*pi))*(cos(za - mu) - cos(za + mu))
% See hand-written notes for details

mu_b = 0.05;                     % half-width of finite support
% mu_b = 0.5*abs(Za - pi/2);

NumZ = numel(z);                % number of z intervals

% colocated actuators and sensors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumAc = 1;
ActPos = Za;

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
Binf = betaU*mu_b*sqrt(1/(2*pi))*(cos(ActPos - mu_b) - cos(ActPos + mu_b));
Bmodel = betaUmodel*sqrt(2/pi)*sin(ActPos);

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

%% a coordinate
a_step = 0.001;
a_BC = 1;
a = -a_BC:a_step:a_BC;
NumA = length(a);


da = zeros(1,NumA);
deltaf = zeros(1,NumA);
ga = zeros(1,NumA);
gamodel = zeros(1,NumA);
deltag = zeros(1,NumA);

Dda = zeros(1,NumA-1);
Ddeltaf = zeros(1,NumA-1);
Dga = zeros(1,NumA-1);
Ddeltag = zeros(1,NumA-1);

for i = 1:NumA    
    
    X = Phi'*Za;
    da(i) = NonlinearFunction(Ainf,a(i),gamma,betaTreal,Phi,X,Zstep);
    
    damodel = NonlinearFunction(Amodel,a(i),gamma,betaT,Phi_model,X,Zstep);
    deltaf(i) = da(i) - damodel;
    
%     da(i) = Ainf*a(i);
%     damodel = Amodel*a(i) + Bmodel*U;
%     deltaf(i) = da(i) - damodel;

    U = -K*a(i);
    ga(i) =  Binf*U;
    gamodel(i) = Bmodel*U;
    deltag(i) = ga(i) - gamodel(i);
end

for i = 1:NumA-1
   Dda(i) = (da(i+1) - da(i))/a_step;
   Ddeltaf(i) = (deltaf(i+1) - deltaf(i))/a_step;
   Dga(i) = (ga(i+1) - ga(i))/a_step;
   Ddeltag(i) = (deltag(i+1) - deltag(i))/a_step;
end

Kf = max(abs(Dda));
Km = max(abs(Dga));
Kdelta = max(abs(Ddeltaf) + abs(Ddeltag));

%%
% figure(1)
% plot(ga);hold on; plot(gamodel, 'r-'); hold on;
% 
% 
% figure()
% plot(a(1:NumA-1), Dda); hold on;
% plot(a(1:NumA-1), abs(Dda), 'r');
% xlim([-a_BC a_BC]);
% legend('f(a)', '|f(a)|');
% xlabel('a','FontSize',30);
% ylabel('f(a)','FontSize',30)
% title('Find $K_f = 0.0074$ from f(a)');
% set(gca,'FontSize',24);
% 
% figure()
% plot(a(1:NumA-1), Dga); hold on;
% plot(a(1:NumA-1), abs(Dga), 'r');
% xlim([-a_BC a_BC]);
% legend('g(a)', '|g(a)|');
% xlabel('a','FontSize',30);
% ylabel('g(a)','FontSize',30)
% title('Find $K_m = 0.0094$ from f(a)');
% set(gca,'FontSize',24);
% 
% figure()
% plot(a(1:NumA-1), Ddeltaf); hold on;
% plot(a(1:NumA-1), abs(Ddeltaf), 'r');
% xlim([-a_BC a_BC]);
% legend('$\delta_f(a)$', '$|\delta_f(a)|$');
% xlabel('$a$','FontSize',30);
% ylabel('$\delta_f(a)$','FontSize',30)
% title('Find $K_{\delta_g} = 6.4837e-04$ from $f(a)$');
% set(gca,'FontSize',24);
% 
% figure()
% plot(a(1:NumA-1), Ddeltag); hold on;
% plot(a(1:NumA-1), abs(Ddeltag), 'r');
% xlim([-a_BC a_BC]);
% legend('$\delta_g(a)$', '$|\delta_g(a)|$');
% xlabel('$a$','FontSize',30);
% ylabel('$\delta_g(a)$','FontSize',30)
% title('Find $K_{\delta_g} = 6.4837e-04$ from $f(a)$');
% set(gca,'FontSize',24);