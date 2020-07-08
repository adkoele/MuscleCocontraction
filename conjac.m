function J = conjac(X, params)
% define constraints jacobian

h = params.h; %time step
N = params.N;

nvarpernode = params.nvarpernode;
nvars = params.nvars;
ncon = params.ncon;
nstates = params.nstates;
ncontrols = params.ncontrols;
ndof = params.ndof;

ix1 = ncontrols+(1:nstates);
iu = 1:ncontrols;
u = X(iu);
ic = 1:nstates;

J = spalloc(ncon,nvars,params.Jnnz);

% Dynamics should match next node till one before last node
for i = 1:N
    if i == N
        ix2 = ncontrols+(1:nstates);
    else
        ix2 = ix1 + nvarpernode;
    end

    x1 = X(ix1);
    x2 = X(ix2);

    %Time delay
    xu1 = x1;
    xu2 = x2;
    ixu1 = ix1;
    ixu2 = ix2;
    
    omega_c = params.c_omega(:,i);
    [u1, du1dx, du1dK, du1dKd] =  findTorque(u,X(end-params.Ks+1:end),xu1(1:ndof*2),params);
    [u2, du2dx, du2dK, du2dKd] =  findTorque(u,X(end-params.Ks+1:end),xu2(1:ndof*2),params);
        
    if strcmp(params.muscletype, 'SEE_linear') || strcmp(params.muscletype, 'SEE_quad')
        du1dx = du1dx';
        du2dx = du2dx'; 
    elseif strcmp(params.muscletype, 'NoAct_linear') || strcmp(params.muscletype, 'NoAct_quad')
        du1dx = [du1dx' zeros(params.nmus,params.nmus)];
        du2dx = [du2dx' zeros(params.nmus,params.nmus)];
    else
        du1dx = [du1dx' zeros(params.nmus,params.nmus*2)];
        du2dx = [du2dx' zeros(params.nmus,params.nmus*2)];
    end
    
    u1 = u1+omega_c;
    u2 = u2+omega_c;

    omega = params.omega(:,i);
    if strcmp(params.method, 'ME')
        [~, dfdx, dfdxdot, dfdu] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);
        dfdu1 = dfdu/2;
        dfdu2 = dfdu/2;
        dfdK = dfdu1/2*du1dK+dfdu2/2*du2dK;
        dfdKd= dfdu1/2*du1dKd+dfdu2/2*du2dKd;


        J(ic,ix1) = dfdx/2 - dfdxdot/h;
        J(ic,ix2) = dfdx/2 + dfdxdot/h;
        if ~isempty(ixu1)
            J(ic,ixu1) = J(ic,ixu1) + dfdu/2*du1dx;
        end
        if ~isempty(ixu2)
            J(ic,ixu2) = J(ic,ixu2) + dfdu/2*du2dx;
        end
    else
        [~, dfdx, dfdxdot, dfdu] = StocDyn(x2,(x2-x1)/h,u2, omega, params);
        dfdu1 = zeros(size(dfdu));
        dfdu2 = dfdu;
        dfdK = dfdu2*du2dK;
        dfdKd= dfdu2*du2dKd;

        J(ic,ix1) = -dfdxdot/h;
        J(ic,ix2) = dfdx + dfdxdot/h;
        if ~isempty(ixu2)
            J(ic,ixu2) = J(ic,ixu2) + dfdu*du2dx;
        end
    end

    J(ic,iu) = (dfdu1+dfdu2);
    J(ic,end-params.Ks+1:end) = [dfdK dfdKd];

    J(end+i-1-params.N,ix1(1:ndof)) = 1;
    
    ix1 = ix1+nvarpernode;
    ic = ic+nstates;
end