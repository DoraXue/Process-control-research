function U = NonlinearController(amodel,LfV,LgV,LwV,theta1)

rho = 0.001;
chi = 2.1;
sai = 1e-3;
e = 1e-15;
if norm(LgV) == 0
    U = zeros(3,1);
else
    na = norm(amodel);

    LfV2 = LfV + rho*na + chi*norm(LwV)*theta1;
    LfV1 = LfV + (LfV2 - LfV)*(na/(na + sai));

    U = -(LfV1 + sqrt(LfV^2 + norm(LgV')^4))/(norm(LgV)^2 + e)*(LgV)';
end
