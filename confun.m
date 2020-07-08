function [c] = confun(X, params)
% define constraints

nstates = params.nstates;
ncontrols = params.ncontrols;
h = params.h;
N = params.N;

ix = ncontrols+(1:nstates);
u = X(1:ncontrols); 
ic = 1:nstates;
        
nvarpernode = params.nvarpernode;
ncon = params.ncon;
ndof = params.ndof;

c = zeros(ncon,1);

% Constraints on dynamics
for i = 1:N
    x1 = X(ix);
    if i == N
        x2 = X(ncontrols+(1:nstates));
    else
        x2 = X(ix+nvarpernode);
    end

    xu1 = x1;
    xu2 = x2;
    
    omega_c = params.c_omega(:,i);
    u1 =  findTorque(u,X(end-params.Ks+1:end),xu1(1:ndof*2),params);
    u2 =  findTorque(u,X(end-params.Ks+1:end),xu2(1:ndof*2),params);
    
    u1 = u1+omega_c;
    u2 = u2+omega_c;
    
    omega = params.omega(:,i);
    if strcmp(params.method, 'ME')
        dyns = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);
    else
        dyns = StocDyn(x2,(x2-x1)/h,u2, omega, params);
    end
    c(ic) = dyns;
    
    c(end+i-1-params.N) = x1(1:ndof)-params.targetangle;
    
    ix = ix+nvarpernode;
    ic = ic+nstates;
end