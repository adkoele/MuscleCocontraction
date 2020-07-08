function params = getParams(params)

params.nvarpernode = params.nstates;
params.nvarSU = params.nvarpernode*params.N;
params.nvars = params.ncontrols+params.nvarSU+params.Ks;
params.ncon = params.nstates*(params.N)+params.N;%+1;