function dobj = objgrad(X, params)

% keyboard;
ix = params.ncontrols+(1:params.nstates);
iact = params.ncontrols+params.ndof*2+(1:params.nmus);
ilce = iact+params.nmus;

dobj = zeros(size(X));
for i = 1:params.N
    [u, dudx, dudK, dudKd] = findTorque(X(1:params.ncontrols),X(end-params.Ks+1:end),X(ix(1:params.ndof*2)),params);
    
    if strcmp(params.muscletype, 'SEE_linear') || strcmp(params.muscletype, 'SEE_quad')
        [u,dudx1,dudu]  = getMusDyns_Fseeonly(X(ix),u,params);
        dudx1 = dudx1./params.muscleparam.fmax';
        dudu = dudu./params.muscleparam.fmax';
        u = u./params.muscleparam.fmax';
    end
    
    if strcmp(params.muscletype, 'SEE_linear') || strcmp(params.muscletype, 'SEE_quad')
        dudx = dudx';
    elseif strcmp(params.muscletype, 'NoAct_linear') || strcmp(params.muscletype, 'NoAct_quad')
        dudx = [dudx' zeros(params.nmus,params.nmus)];
    else
        dudx = [dudx' zeros(params.nmus,params.nmus*2)];
    end

    % Muscles
    if strcmp(params.muscletype, 'SEE_linear') || strcmp(params.muscletype, 'SEE_quad')
        dobj(1:params.ncontrols) = dobj(1:params.ncontrols) + dudu*u;
        dobj(ix) = dobj(ix) + dudx1'*u + dudx'*dudu*u;%
        dobj(end-params.Ks+1:end) = dobj(end-params.Ks+1:end)+[dudK*dudu*u;dudKd*dudu*u];
    else
        dobj(ix) = dobj(ix) + dudx'*u;
        
        dobj(end-params.Ks+1:end) = dobj(end-params.Ks+1:end)+[dudK*u;dudKd*u];
        dobj(1:params.ncontrols) = dobj(1:params.ncontrols)+u;
    end
    ix = ix+params.nvarpernode;
    ilce = ilce+params.nvarpernode;
    iact = iact+params.nvarpernode;
end
dobj = dobj/params.N/params.stdev*100;