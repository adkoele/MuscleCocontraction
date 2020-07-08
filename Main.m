rng('shuffle')

clear all
close all
clc

%Main program for pendulum swingup

% Declare parameters
params.m = 2;
params.l = 0.6;
params.g = 9.81;
params.T = 3; %Duration of motion
params.N = 3000;
params.h = params.T/(params.N-1);

params.ndof = 1;
params.nmus = 2;
params.ncontrols= params.nmus;
params.Ks = params.ncontrols*2;
params.targetangle = pi/2;

params.muscleparam = Muscleparams();
params.solver = 'IPOPT';
params.method = 'BE';
params.epsilon = 1e-2;

params.muscletype = 'quad'; %options: 'quad', 'linear', 'SEE_linear', 'SEE_quad', 'NoAct_quad', 'NoAct_linear'

if contains(params.muscletype, 'SEE')
    params.nstates = params.ndof*2;
elseif contains(params.muscletype, 'NoAct')
    params.nstates = params.ndof*2+params.nmus;
else
    params.nstates = params.ndof*2+params.nmus*2;
end

params.xneutral = findNeutralstate(params);

lengthi = 1; %noise repetitions
theta_bounds = [5 10 15 20]/180*pi;
stdevs = [0 2 5 10 15 20]; 

params.controlnoise = 0;
params.inputnoise = 0;

for ii = 1:lengthi
    randvals = rand(params.ndof*2,params.N);
    for k = 1:length(theta_bounds)
        params.theta_bound = theta_bounds(k);
        params.stdev = 1;
        
        %First result witout noise for desired trajectory
        params = getParams(params);
        params.omega = zeros(params.ndof*2,params.N); %Added noise
        params.c_omega = params.omega;
        [X0, L, U] = getIniConBound(params,1);
        [~, params] = getConjacstructure(L, U, params);
        
        %Derivative check
%         [grad,grad_num,cjac, cjac_num] = checkDerivatives(X0+randn(size(X0))*0.1, params);%
%         keyboard;

        params.warmstart = 0;
        result1 = Optimize(X0, L, U, params);
        result(1) = result1;
        allresult(1,k).result(ii) = result1;
        clear mex

        for i = 2:length(stdevs)
            params.warmstart = 1;

            % Now add noise and increase time
            stdev = stdevs(i)/sqrt(params.h);
            params.stdev = stdev;
            
            if ~params.controlnoise
                params.omega = -stdev/2 + stdev*randvals;
                params.c_omega = zeros(size(randvals));
            else
                params.c_omega = -stdev/2 + stdev*randvals;
                params.omega = zeros(size(randvals));
            end

            [X0, L, U] = getIniConBound(params,1,result(i-1));%result1);
            [~, params] = getConjacstructure(L, U, params);

            if strcmp(params.solver,'IPOPT')
                params.zl = result(i-1).zl;
                params.zu = result(i-1).zu;
                params.lambda = result(i-1).lambda;
            end

            result2 = Optimize(X0, L, U, params);
            result(i) = result2;
            allresult(i,k).result(ii) = result2;
        end
    end
end
