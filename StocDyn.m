function [f, dfdx, dfdxdot, dfdu] = StocDyn(x, xdot, u, omega, params)

m = params.m;
l = params.l;
g = params.g;
nmus = params.nmus;
ndof = params.ndof;

theta = x(1);
dtheta = x(2);
b = 0;

%Run muscle dynamics
if strcmp(params.muscletype, 'linear')
    [Fsee,dFseedx,fmus,dfmusdx,dfmusdxdot,dfmusdu] = getMusDyns_linear(x,xdot,u,params);
elseif strcmp(params.muscletype, 'SEE_linear') || strcmp(params.muscletype, 'SEE_quad')
    [Fsee,dFseedx,dFseedlce]  = getMusDyns_Fseeonly(x,u,params);
elseif strcmp(params.muscletype, 'NoAct_linear') || strcmp(params.muscletype, 'NoAct_quad')
    [Fsee,dFseedx,fmus,dfmusdx,dfmusdxdot,dfmusdu] = getMusDyns_NoAct(x,xdot,u,params);
else
    [Fsee,dFseedx,fmus,dfmusdx,dfmusdxdot,dfmusdu] = getMusDyns(x,xdot,u,params);
end

% Find torque
d = params.muscleparam.d;
T = d*Fsee;
dTdx = d*dFseedx;

J = m*l^2;
G = m*g*l*cos(theta);

a_x = [dtheta; -G/J];
Bu = [0; (T-b*dtheta)/J];
C = zeros(params.ndof*2);
C(2,2) = 1;

f = a_x + Bu + C*omega - xdot(1:params.ndof*2);

if strcmp(params.muscletype, 'SEE_linear') || strcmp(params.muscletype, 'SEE_quad')
    dfdx = [0 1; m*g*l*sin(theta)/J -b/J];
    dfdx(ndof*2,:) = dfdx(ndof*2,:)+dTdx/J;
    dfdxdot = -1*eye(params.ndof*2);
    dfdu = [0 0; d*dFseedlce/J];
elseif strcmp(params.muscletype, 'NoAct_linear') || strcmp(params.muscletype, 'NoAct_quad')
    f = [f;fmus];
    dfdx = [[0 1; m*g*l*sin(theta)/J -b/J] zeros(params.ndof*2,params.nmus);dfmusdx];
    dfdx(ndof*2,:) = dfdx(ndof*2,:)+dTdx/J;
    dfdu = [zeros(ndof*2,params.nmus);dfmusdu];
    dfdxdot = [-1*eye(params.ndof*2) zeros(params.ndof*2,params.nmus);dfmusdxdot];
else
    f = [f;fmus];
    dfdx = [[0 1; m*g*l*sin(theta)/J -b/J] zeros(params.ndof*2,params.nmus*2);dfmusdx];
    dfdx(ndof*2,:) = dfdx(ndof*2,:)+dTdx/J;
    dfdu = [zeros(ndof*2,params.nmus);dfmusdu];
    dfdxdot = [-1*eye(params.ndof*2) zeros(params.ndof*2,params.nmus*2);dfmusdxdot];
end

%scaling
factor = 1e3;
f = f/factor;
dfdx = dfdx/factor;
dfdu = dfdu/factor;
dfdxdot = dfdxdot/factor;