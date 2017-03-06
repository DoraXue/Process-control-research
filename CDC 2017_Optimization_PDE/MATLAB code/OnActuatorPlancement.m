clear;clc;

% construct z axis
BC_ini = 0;                     % boundary condition 1
BC_fin = pi;                    % boundary condition 2
Zstep = pi/100;                 % step of z
z = BC_ini: Zstep: BC_fin;
NumZ = numel(z);                % number of z intervals

NumEv= 3;
NumAc = 3;

% other constant parameters
A = [1 0 0;
     0 -2 0;
     0 0 -5];
deltaA = [0.02 0 0;
          0 0.02 0;
          0 0 0.02];
Amodel = A + deltaA;

deltaB = 0.5*deltaA;
P = eye(3);
poles = [-100 -100 -100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Case1: 1 degree of free location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ActPos = zeros(1,NumAc);
ActPos(1,1) = pi/3;
ActPos(1,2) = 2*pi/3;

Thresh = zeros(1,NumZ);
ControllableRegion = zeros(1,NumZ);
CH = zeros(1,NumZ);
alpha = zeros(1,NumZ);

for loc = 1 : NumZ
    Act3 = (loc-1)*Zstep;
    ActPos(1,3) = Act3;
    SenPos = ActPos;                % colocated sensors and actuators
    % C matrix
    Q = zeros(NumEv,NumAc);
    for n = 1: NumEv
        for i = 1:NumAc
        Q(n,i) = sqrt(2/pi)*sin(n*SenPos(1,i)); 
        end
    end
    C = Q';
    
    % B matrix
    B = zeros(NumEv,NumAc);
    for n = 1: NumEv
        for i = 1:NumAc
            B(n,i) = 2*sqrt(2/pi)*sin(n*ActPos(1,i)); 
        end
    end
    Bmodel = B + deltaB;
    
    % check controllability
    Co = ctrb(A,B);
    if rank(Co) == 3
        ControllableRegion(1,loc) = 1;      % controllable
    else
        ControllableRegion(1,loc) = 0;      % uncontrollable
    end
    
    % check C's invertability
    detC = det(C);
    if detC == 0
        CH(1,loc) = 0;
    else
        CH(1,loc) = 1;
    end
    
    % design controller
    K = place(Amodel,Bmodel,poles);
    kk = C^(-1)*K;
    
    % calculate alpha
    Q = -((Amodel - Bmodel*K)'*P + P*(Amodel - Bmodel*K));
    alpha(1,loc) = min(eig(Q));
    
    % calculate
    Thresh(1,loc) = (alpha(1,loc)/2 - 2*norm(P*deltaA) - 2*norm(P*deltaB)*norm(K))/(2*norm(P*(Bmodel*norm(K) + deltaB*norm(K))));

end
% plot the maginitude of threshold and actuator locations
figure()
plot(z,Thresh)
axis([0 pi -0.1 0.4])
xlabel('Z','FontSize',30);
ylabel('Error threshold coefficient','FontSize',30)
set(gca,'FontSize',24);

% plot controllable regions
figure()
plot(z,ControllableRegion)
axis([0 pi -0.5 1.5])
xlabel('Z','FontSize',30);
ylabel('Controllable region','FontSize',30)
set(gca,'FontSize',24);

% plot alpha
figure()
plot(z,alpha)
xlim = ([0 pi]);
% axis([0 pi 100 300])
xlabel('Z','FontSize',30);
ylabel('alpha','FontSize',30)
set(gca,'FontSize',24);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Case 2: 2 degrees of free locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ActPos = zeros(1,NumAc);
% ActPos(1,1) = pi/4;
% 
% Thresh = zeros(NumZ,NumZ);
% ControllableRegion = zeros(NumZ,NumZ);
% CH = zeros(NumZ,NumZ);
% alpha = zeros(NumZ,NumZ);
% 
% for z2 = 1 : NumZ
%     Act2 = (z2-1)*Zstep;
%     ActPos(1,2) = Act2;
%    for z3 = 1 : NumZ
%        Act3 = (z3-1)*Zstep;
%        ActPos(1,3) = Act3;
%        SenPos = ActPos;
%        
%        % C matrix
%     Q = zeros(NumEv,NumAc);
%     for n = 1: NumEv
%         for i = 1:NumAc
%         Q(n,i) = sqrt(2/pi)*sin(n*SenPos(1,i)); 
%         end
%     end
%     C = Q';
%     
%     % B matrix
%     B = zeros(NumEv,NumAc);
%     for n = 1: NumEv
%         for i = 1:NumAc
%             B(n,i) = 2*sqrt(2/pi)*sin(n*ActPos(1,i)); 
%         end
%     end
%     Bmodel = B + deltaB;
%     
%     % check controllability
%     Co = ctrb(A,B);
%     if rank(Co) == 3
%         ControllableRegion(z2,z3) = 1;  % controllable
%     else
%         ControllableRegion(z2,z3) = 0;  % uncontrollalbe
%     end
%     
%     % check C's invertability
%     detC = det(C);
%     if detC == 0
%         CH(z2,z3) = 1;                  % invertable
%     else
%         CH(z2,z3) = 0;                  % non-invertable
%     end
%     
%     % design controller
%     K = place(Amodel,Bmodel,poles);
%     kk = C^(-1)*K;
%     
%     % calculate alpha
%     Q = -((Amodel - Bmodel*K)'*P + P*(Amodel - Bmodel*K));
%     alpha(z2,z3) = min(eig(Q));
%     
%     % calculate
%     Thresh(z2,z3) = (alpha(z2,z3)/2 - 2*norm(P*deltaA) - 2*norm(P*deltaB)*norm(K))/(2*norm(P*(Bmodel*norm(K) + deltaB*norm(K))));
%    end
% end
% 
% % plot 3D profile of threshold
% figure()
% mesh(z,z,Thresh);
% axis([0 pi 0 pi -0.1 0.5])
% xlabel('Z1','FontSize',30);
% ylabel('Z2','FontSize',30)
% zlabel('Error threshold coefficient','FontSize',30);
% set(gca,'FontSize',24);
% 
% % contour plot of threshold
% figure()
% contour(z,z,Thresh)
% [C_label,h_label] = contour(z,z,Thresh);
% clabel(C_label,h_label);
% xlabel('Z1','FontSize',30);
% ylabel('Z2','FontSize',30)
% zlabel('Error threshold coefficient','FontSize',30);
% set(gca,'FontSize',24);
% 
% % plot controllable region
% figure()
% mesh(z,z,ControllableRegion);
% % axis([0 pi 0 pi 0.08 0.2])
% xlabel('Z1','FontSize',30);
% ylabel('Z2','FontSize',30)
% zlabel('Controllable region','FontSize',30);
% set(gca,'FontSize',24);
% 
% % plot alpha
% figure()
% mesh(z,z,alpha);
% % axis([0 pi 0 pi 0.08 0.2])
% xlabel('Z1','FontSize',30);
% ylabel('Z2','FontSize',30)
% zlabel('Alpha','FontSize',30);
% set(gca,'FontSize',24);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Case 3: fix all sensors; 1 degree of free actuator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ActPos = zeros(1,NumAc);
% ActPos(1,1) = pi/3;
% ActPos(1,2) = 2*pi/3;
% 
% SenPos = zeros(1,NumAc);
% SenPos(1,1) = pi/4;
% SenPos(1,2) = 2*pi/4;
% SenPos(1,3) = 3*pi/4;
% 
% Thresh = zeros(1,NumZ);
% ControllableRegion = zeros(1,NumZ);
% CH = zeros(1,NumZ);
% alpha = zeros(1,NumZ);
% 
% for loc = 1 : NumZ
%     Act3 = (loc-1)*Zstep;
%     ActPos(1,3) = Act3;
%  
%     % C matrix
%     Q = zeros(NumEv,NumAc);
%     for n = 1: NumEv
%         for i = 1:NumAc
%         Q(n,i) = sqrt(2/pi)*sin(n*SenPos(1,i)); 
%         end
%     end
%     C = Q';
%     
%     % B matrix
%     B = zeros(NumEv,NumAc);
%     for n = 1: NumEv
%         for i = 1:NumAc
%             B(n,i) = 2*sqrt(2/pi)*sin(n*ActPos(1,i)); 
%         end
%     end
%     Bmodel = B + deltaB;
%     
%     % check controllability
%     Co = ctrb(A,B);
%     if rank(Co) == 3
%         ControllableRegion(1,loc) = 1;      % controllable
%     else
%         ControllableRegion(1,loc) = 0;      % uncontrollable
%     end
%     
%     % check C's invertability
%     detC = det(C);
%     if detC == 0
%         CH(1,loc) = 0;
%     else
%         CH(1,loc) = 1;
%     end
%     
%     % design controller
%     K = place(Amodel,Bmodel,poles);
%     kk = C^(-1)*K;
%     
%     % calculate alpha
%     Q = -((Amodel - Bmodel*K)'*P + P*(Amodel - Bmodel*K));
%     alpha(1,loc) = min(eig(Q));
%     
%     % calculate
%     Thresh(1,loc) = (alpha(1,loc)/2 - 2*norm(P*deltaA) - 2*norm(P*deltaB)*norm(K))/(2*norm(P*(Bmodel*norm(K) + deltaB*norm(K))));
% 
% end
% % plot the maginitude of threshold and actuator locations
% figure()
% plot(z,Thresh)
% axis([0 pi -0.05 0.3])
% xlabel('Z','FontSize',30);
% ylabel('Error threshold coefficient','FontSize',30)
% set(gca,'FontSize',24);

% % plot controllable regions
% figure()
% plot(z,ControllableRegion)
% axis([0 pi -0.5 1.5])
% xlabel('Z','FontSize',30);
% ylabel('Controllable region','FontSize',30)
% set(gca,'FontSize',24);
% 
% % plot alpha
% figure()
% plot(z,alpha)
% xlim = ([0 pi]);
% % axis([0 pi 100 300])
% xlabel('Z','FontSize',30);
% ylabel('alpha','FontSize',30)
% set(gca,'FontSize',24);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Case 4: 2 degrees of free locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ActPos = zeros(3,NumAc);
% ActPos(1,1) = pi/2;
% 
% SenPos = zeros(1,NumAc);
% SenPos(1,1) = pi/4;
% SenPos(1,2) = 2*pi/4;
% SenPos(1,3) = 3*pi/4;
% 
% Thresh = zeros(NumZ,NumZ);
% ControllableRegion = zeros(NumZ,NumZ);
% CH = zeros(NumZ,NumZ);
% alpha = zeros(NumZ,NumZ);
% 
% for z2 = 1 : NumZ
%     Act2 = (z2-1)*Zstep;
%     ActPos(1,2) = Act2;
%    for z3 = 1 : NumZ
%        Act3 = (z3-1)*Zstep;
%        ActPos(1,3) = Act3;
%        
%        % C matrix
%     Q = zeros(NumEv,NumAc);
%     for n = 1: NumEv
%         for i = 1:NumAc
%         Q(n,i) = sqrt(2/pi)*sin(n*SenPos(1,i)); 
%         end
%     end
%     C = Q';
%     
%     % B matrix
%     B = zeros(NumEv,NumAc);
%     for n = 1: NumEv
%         for i = 1:NumAc
%             B(n,i) = 2*sqrt(2/pi)*sin(n*ActPos(1,i)); 
%         end
%     end
%     Bmodel = B;
%     
%     % check controllability
%     Co = ctrb(A,B);
%     if rank(Co) == 3
%         ControllableRegion(z2,z3) = 1;  % controllable
%     else
%         ControllableRegion(z2,z3) = 0;  % uncontrollalbe
%     end
%     
%     % check C's invertability
%     detC = det(C);
%     if detC == 0
%         CH(z2,z3) = 1;                  % invertable
%     else
%         CH(z2,z3) = 0;                  % non-invertable
%     end
%     
%     % design controller
%     K = place(Amodel,Bmodel,poles);
%     kk = C^(-1)*K;
%     
%     % calculate alpha
%     Q = -((Amodel - Bmodel*K)'*P + P*(Amodel - Bmodel*K));
%     alpha(z2,z3) = min(eig(Q));
%     
%     % calculate
%     Thresh(z2,z3) = (alpha(z2,z3)/2 - 2*norm(P*deltaA) - 2*norm(P*deltaB)*norm(K))/(2*norm(P*(Bmodel*norm(K) + deltaB*norm(K))));
%    end
% end
% 
% % plot 3D profile of threshold
% figure()
% mesh(z,z,Thresh);
% axis([0 pi 0 pi -0.1 0.5])
% xlabel('Z1','FontSize',30);
% ylabel('Z2','FontSize',30)
% zlabel('Error threshold coefficient','FontSize',30);
% set(gca,'FontSize',24);
% 
% % contour plot of threshold
% figure()
% contour(z,z,Thresh)
% [C_label,h_label] = contour(z,z,Thresh);
% clabel(C_label,h_label);
% xlabel('Z1','FontSize',30);
% ylabel('Z2','FontSize',30)
% zlabel('Error threshold coefficient','FontSize',30);
% set(gca,'FontSize',24);
% 
% % plot controllable region
% figure()
% mesh(z,z,ControllableRegion);
% % axis([0 pi 0 pi 0.08 0.2])
% xlabel('Z1','FontSize',30);
% ylabel('Z2','FontSize',30)
% zlabel('Controllable region','FontSize',30);
% set(gca,'FontSize',24);
% 
% % plot alpha
% figure()
% mesh(z,z,alpha);
% % axis([0 pi 0 pi 0.08 0.2])
% xlabel('Z1','FontSize',30);
% ylabel('Z2','FontSize',30)
% zlabel('Alpha','FontSize',30);
% set(gca,'FontSize',24);
% 

