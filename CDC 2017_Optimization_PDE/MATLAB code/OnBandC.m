clear;clc;

% construct z axis
BC_ini = 0;                     % boundary condition 1
BC_fin = pi;                    % boundary condition 2
Zstep = pi/100;                 % step of z
z = BC_ini: Zstep: BC_fin;
NumZ = numel(z);                % number of z intervals

ValidLoc = zeros(1,NumZ);

for z1 = 1: NumZ
   za = (z1-1)*Zstep;
   B = sin(za);
   C = inv(B);
   if C <=1 && C >=0
       
       ValidLoc(1,z1) = 1;
   end
end
plot(z,ValidLoc.'.')
axis([0 pi -0.1 1.1])
xlabel('Z','FontSize',30);
ylabel('Valid actuator location','FontSize',30)
set(gca,'FontSize',24);

% ValidLoc = zeros(NumZ,NumZ);

% for z1 = 1:NumZ
%    za1 = (z1-1)*Zstep;
%     for z2 = 1:NumZ
%         za2 = (z2-1)*Zstep;
%         
%         B = [2*sqrt(2/pi)*sin(za1) 2*sqrt(2/pi)*sin(za2);
%              2*sqrt(2/pi)*sin(2*za1) 2*sqrt(2/pi)*sin(2*za2)];
%         C = inv(B);
%         zs1 = asin(C(1,1)/sqrt(2/pi));
%         zs2 = asin(C(2,1)/sqrt(2/pi));
%         
%         if abs(sqrt(2/pi)*sin(2*zs1) - C(1,2)) < 1e-2
%            if abs(sqrt(2/pi)*sin(2*zs2) - C(2,2)) < 1e-2
%                ValidLoc(z1,z2) = 1;
%            end
%         end
%     end
% end
% 
% for z1 = 1:NumZ
%    za1 = (z1-1)*Zstep;
%     for z2 = 1:NumZ
%         za2 = (z2-1)*Zstep;
%         
%         B = [sin(za1) sin(za2);
%              sin(2*za1) sin(2*za2)];
%         C = inv(B);
%         zs1 = asin(C(1,1));
%         zs2 = asin(C(2,1));
%         
%         if abs(sin(2*zs1) - C(1,2)) < 1e-2
%            if abs(sin(2*zs2) - C(2,2)) < 1e-2
%                ValidLoc(z1,z2) = 1;
%            end
%         end
%     end
% end
% 
% figure()
% mesh(z,z,ValidLoc);
% axis([0 pi 0 pi -0.1 0.5])
% xlabel('Z1','FontSize',30);
% ylabel('Z2','FontSize',30)
% zlabel('Valid actuator locations','FontSize',30);
% set(gca,'FontSize',24);