% Explored why the stability region has no dependency on Za
% It was because the ratio of B/Bmodel is constant, a function of mu

clear;clc;

BC_ini = 0;                     % boundary condition 1
BC_fin = pi;                    % boundary condition 2
Zstep = pi/100;                 % step of z
z = BC_ini: Zstep: BC_fin;
NumZ = numel(z);    

B = zeros(1,NumZ);
Bhat = zeros(1,NumZ);
ratio = zeros(1,NumZ);

mu = 0.1;

for i = 1:NumZ
    za = z(i);
    B(i) = mu*sqrt(1/(2*pi))*(cos(za - mu) - cos(za + mu));
    Bhat(i) = sqrt(2/pi)*sin(za);
    ratio(i) = B(i)/Bhat(i);
end
%%
figure()
plot(z, B); hold on;
plot(z, Bhat, 'r')
xlim([0 pi]);
legend('$B(z_a)$', '$\hat{B}(z_a)$');
xlabel('$z_a$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',24);

figure()
plot(z, ratio)
xlim([0 pi]);
xlabel('$z_a$','Interpreter','latex','FontSize',30);
ylabel('$B(z_a)/\hat{B}(z_a)$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',24);
