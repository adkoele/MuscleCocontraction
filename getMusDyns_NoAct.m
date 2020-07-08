function [Fsee,dFseedx,f,dfdx,dfdxdot,dfdu, stiffness] = getMusDyns_NoAct(x,xdot,u,params)

%muscle parameters
A = params.muscleparam.A;
vmax = params.muscleparam.vmax;
w = params.muscleparam.w;
d = params.muscleparam.d';
fmax = params.muscleparam.fmax';
gmax = params.muscleparam.gmax;
lceopt = params.muscleparam.lceopt';
lsees = params.muscleparam.lsees';
l0 = params.muscleparam.l0';
ndof = params.ndof;
nmus = params.nmus;
nstates = params.nstates;
ksee = params.muscleparam.ksee';

%Initialize f dfdx dfdxdot dfdu
f = zeros(nmus,1);
dfdx = zeros(nmus,nstates);
dfdxdot = zeros(nmus,nstates);

%saturate u
a = 1/2*(u+sqrt(u.^2+10^-4));
dadu = 1/2*(1+u./sqrt(u.^2+10^-4));

%states
angle = (x(1:ndof)-pi/2);
lce = x(2*ndof+(1:nmus));
lcedot = xdot(2*ndof+(1:nmus));

%Force in contractile element
flce = exp(-((lce-1)/w).^2);
dflcedlce = -2*(lce-1)/w^2.*flce;
glce = zeros(nmus,1);
dglcedlcedot = zeros(nmus,1);
for i = 1:nmus
    if lcedot(i) < 0
        glce(i) = (vmax+lcedot(i))/(vmax-lcedot(i)/A);
        dglcedlcedot(i) = 1/(vmax-lcedot(i)/A)+(vmax+lcedot(i))/((vmax-lcedot(i)/A)^2*A);
    else
        c3 = vmax*A*(gmax-1)/(A+1);
        glce(i) = (gmax*lcedot(i)+c3)/(lcedot(i)+c3);
        dglcedlcedot(i) = gmax/(lcedot(i)+c3)-(gmax*lcedot(i)+c3)/(lcedot(i)+c3)^2;
    end
end
Flce = a.*fmax.*flce.*glce+lcedot.*fmax/1000; % damping force
dFlcedx = [zeros(nmus,ndof*2) diag(a.*fmax.*glce.*dflcedlce)];
dFlcedxdot = [zeros(nmus,ndof*2) diag(a.*fmax.*flce.*dglcedlcedot+fmax/1000)];

% %Force in the parallel elastic element
% Fpee = 0.1*fmax.*lceopt.*(lce-lpees)+(lce-lpees>0).*kpee.*fmax.*lceopt.*(lce-lpees).^2;
% dFpeedx = [zeros(nmus,ndof*2+nmus) diag(0.1*fmax.*lceopt+(lce-lpees>0)*2.*kpee.*fmax.*lceopt.*(lce-lpees))]; 

%Force in the series elastic element, using total muscle length
L_m = l0-angle*d;
dL_mdangle = -d;
Fsee = 0.1*fmax.*(L_m-lce.*lceopt-lsees)+(L_m-lce.*lceopt-lsees > 0).*ksee.*fmax.*(L_m-lce.*lceopt-lsees).^2;
dFseedx = [0.1*fmax.*dL_mdangle+(L_m-lce.*lceopt-lsees > 0)*2.*ksee.*fmax.*(L_m-lce.*lceopt-lsees).*dL_mdangle zeros(nmus,ndof) ...
    diag(-0.1*fmax.*lceopt-(L_m-lce.*lceopt-lsees>0)*2.*ksee.*fmax.*(L_m-lce.*lceopt-lsees).*lceopt)];

% %Imbalance on lce
% f(nmus+1:nmus*2) = (Fsee-Fpee-Flce)./fmax;
% dfdx(nmus+1:nmus*2,:) = bsxfun(@rdivide,(dFseedx-dFpeedx-dFlcedx),fmax);
% dfdxdot(nmus+1:nmus*2,:) = bsxfun(@rdivide,-dFlcedxdot,fmax);
%Imbalance on lce
f(1:nmus) = (Fsee-Flce)./fmax;
dfdu = -diag(flce.*glce).*dadu;
dfdx(1:nmus,:) = bsxfun(@rdivide,(dFseedx-dFlcedx),fmax);
dfdxdot(1:nmus,:) = bsxfun(@rdivide,-dFlcedxdot,fmax);
stiffness = (L_m-lce.*lceopt-lsees > 0).*ksee.*fmax.*(L_m-lce.*lceopt-lsees).*2;
% %% HACK return Flce+Fpee
% Fsee = Flce;