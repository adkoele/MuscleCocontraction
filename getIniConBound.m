function [X, L, U] = getIniConBound(params, nou0, result)

if nargin == 1
    nou0 = 0;
end

% Bounds on state

if strcmp(params.muscletype, 'SEE_linear') || strcmp(params.muscletype, 'SEE_quad')
    Ux = [6*pi; 1000];
    Lx = [-5*pi;-1000];
elseif strcmp(params.muscletype, 'NoAct_linear') || strcmp(params.muscletype, 'NoAct_quad') 
    Ux = [6*pi; 1000;1.9;1.9];
    Lx = [-5*pi;-1000;0.1;0.1];
else
    Ux = [6*pi; 1000;40;40;1.9;1.9];
    Lx = [-5*pi;-1000;0;0;0.1;0.1];
end
X0x = params.xneutral;

% Bounds on controls, with or without co-contraction 
if nou0 ~= 1
    Lu = nou0+zeros(params.ncontrols,1);
    Uu = nou0+repmat(1e-4,params.ncontrols,1);
    X0u = 1e-4+nou0+zeros(params.ncontrols,1);
elseif strcmp(params.muscletype, 'SEE_linear') || strcmp(params.muscletype, 'SEE_quad')
    Lu = [0.5; 0.5];
    Uu = [1.1; 1.1];
    X0u = [0.999; 0.999];
else
    Lu = zeros(params.ncontrols,1);
    Uu = repmat(40,params.ncontrols,1);
    X0u = 1e-4+zeros(params.ncontrols,1);
end


L = [Lu; repmat(Lx,params.N,1)];
U = [Uu; repmat(Ux,params.N,1)];
X0 = [X0u;repmat(X0x,params.N,1)];

if params.omega == 0
    L = [L;zeros(params.Ks,1)]; %K and Kd
    U = [U;zeros(params.Ks,1)];
    X = [X0;zeros(params.Ks,1)];
else
    X = [X0;zeros(params.Ks,1)];
    L = [L;repmat([-100;-100],params.Ks/2,1)];
    U = [U;repmat([100;100],params.Ks/2,1)];
end

if nargin == 3
    if result.params.N == params.N
        X = result.X;
    else
        X(1:params.ncontrols) = result.X(1:params.ncontrols);
        X(end-params.Ks+1:end) = result.X(end-params.Ks+1:end);
    end
end
