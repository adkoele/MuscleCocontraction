function [obj] = objfun(X, params)

itheta = params.ncontrols+params.ndof;
iact = params.ncontrols+params.ndof*2+(1:params.nmus);
ilce = iact+params.nmus;

obj = 0;
for i = 1:params.N
    %Time delay
    u = findTorque(X(1:params.ncontrols),X(end-params.Ks+1:end),X((itheta:itheta+1)),params);
        
    if strcmp(params.muscletype, 'SEE_linear') || strcmp(params.muscletype, 'SEE_quad')
        u  = getMusDyns_Fseeonly(X(itheta:itheta+1),u,params)./params.muscleparam.fmax';
    end
    obj_eff = 1/2*sum(u.^2);
    obj = obj+obj_eff;

    itheta = itheta+params.nvarpernode;
    ilce = ilce+params.nvarpernode;
    iact = iact+params.nvarpernode;
end
%The mean activation
obj = obj/params.N/params.stdev*100;
pause(0)