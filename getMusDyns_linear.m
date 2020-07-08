function [Fsee,dFseedx,f,dfdx,dfdxdot,dfdu, stiffness] = getMusDyns_linear(x,xdot,u,params)

%muscle parameters
A = params.muscleparam.A;
vmax = params.muscleparam.vmax;
w = params.muscleparam.w;
Ts = params.muscleparam.Ts;
d = params.muscleparam.d';
fmax = params.muscleparam.fmax';
gmax = params.muscleparam.gmax;
lceopt = params.muscleparam.lceopt';
lsees = params.muscleparam.lsees';
lpees = params.muscleparam.lpees';
l0 = params.muscleparam.l0';
ndof = params.ndof;
nmus = params.nmus;
nstates = params.nstates;
kpee = params.muscleparam.kpee';
ksee = params.muscleparam.ksee';
epsilon = params.epsilon;

%Initialize f dfdx dfdxdot dfdu
f = zeros(nmus*2,1);
dfdx = zeros(nmus*2,nstates);
dfdxdot = zeros(nmus*2,nstates);
dfdu = zeros(nmus*2,nmus);

%saturate u
usat2 = 1/2*(u+sqrt(u.^2+10^-4));
dusat2du = 1/2*(1+u./sqrt(u.^2+10^-4));

%states
angle = (x(ndof)-pi/2);
a = x(2*ndof+(1:nmus));
lce = x(2*ndof+nmus+(1:nmus));
lcedot = xdot(2*ndof+nmus+(1:nmus));
adot = xdot(2*ndof+(1:nmus));

%Imbalance on activation dynamics
c2 = 1/Ts(2);
c1 = 1/Ts(1)-c2;
f(1:nmus) = (usat2-a).*(c1*usat2+c2)-adot;
dfdx(1:nmus,2*ndof+(1:nmus)) = diag(-(c1*usat2+c2)); %params.ndof*2+(1:nmus)
dfdu(1:nmus,:) = (diag(c1*usat2+c2)+diag((usat2-a)*c1)).*dusat2du;
dfdxdot(1:nmus,2*ndof+(1:nmus)) = -eye(nmus);

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
dFlcedx = [zeros(nmus,ndof*2) diag(fmax.*flce.*glce) diag(a.*fmax.*glce.*dflcedlce)];
dFlcedxdot = [zeros(nmus,ndof*2+nmus) diag(a.*fmax.*flce.*dglcedlcedot+fmax/1000)];

% %Force in the parallel elastic element
% Fpee = 0.1*fmax.*lceopt.*(lce-lpees)+(lce-lpees>0).*kpee.*fmax.*lceopt.*(lce-lpees).^2;
% dFpeedx = [zeros(nmus,ndof*2+nmus) diag(0.1*fmax.*lceopt+(lce-lpees>0)*2.*kpee.*fmax.*lceopt.*(lce-lpees))]; 

%Force in the series elastic element, using total muscle length
% %% HACK return Flce+Fpee
% Fsee = Flce;
L_m = l0-angle*d;
dL_mdangle = -d;
ksee_lin = ksee*0.005;
Fsee1 = ksee_lin.*fmax.*(L_m-lce.*lceopt - lsees);
dFsee1dx = [ksee_lin.*fmax.*dL_mdangle zeros(nmus, nmus+ndof) diag(-ksee_lin.*fmax.*lceopt)];
Fsee = 0.1*fmax.*(L_m-lce.*lceopt-lsees) +  1/2*(Fsee1+sqrt(Fsee1.^2+epsilon^2));
dFseedFsee1 = 1/2*(1+Fsee1./sqrt(Fsee1.^2+epsilon^2));
dFseedx =   [0.1*fmax.*dL_mdangle zeros(nmus, nmus+ndof) diag(-0.1*fmax.*lceopt)] + bsxfun(@times,dFseedFsee1,dFsee1dx);
% 
% %Imbalance on lce
% f(nmus+1:nmus*2) = (Fsee-Fpee-Flce)./fmax;
% dfdx(nmus+1:nmus*2,:) = bsxfun(@rdivide,(dFseedx-dFpeedx-dFlcedx),fmax);
% dfdxdot(nmus+1:nmus*2,:) = bsxfun(@rdivide,-dFlcedxdot,fmax);

%Imbalance on lce
f(nmus+1:nmus*2) = (Fsee-Flce)./fmax;
dfdx(nmus+1:nmus*2,:) = bsxfun(@rdivide,(dFseedx-dFlcedx),fmax);
dfdxdot(nmus+1:nmus*2,:) = bsxfun(@rdivide,-dFlcedxdot,fmax);
stiffness = (L_m-lce.*lceopt-lsees > 0).*fmax.*ksee_lin;
% %% HACK return Flce+Fpee
% Fsee = Flce;