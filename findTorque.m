function [u, dudx, dudK, dudKd] = findTorque(uin,Ks,x, params)

Kfactor = 100;
Kdfactor = 10;
if params.ncontrols == 1
    K = sign(params.muscleparam.d')*Ks(1);
    Kd= sign(params.muscleparam.d')*Ks(2);
else
    K = Ks(1:2)*Kfactor;
    Kd= Ks(3:4)*Kdfactor;
end

u = uin+[K Kd]*(x-[params.targetangle;0]);

dudx = [K'; Kd'];
dudK = diag((x(1)-params.targetangle)*[1;1])*Kfactor;
dudKd= diag(x(2)*[1;1])*Kdfactor;
