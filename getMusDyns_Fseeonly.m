function [Fsee,dFseedx,dFseedlce, stiffness]  = getMusDyns_Fseeonly(x,lce,params)

%muscle parameters
lsees = params.muscleparam.lsees';
l0 = params.muscleparam.l0';
lceopt = params.muscleparam.lceopt';
ksee = params.muscleparam.ksee';
fmax = params.muscleparam.fmax';
d = params.muscleparam.d';
epsilon = params.epsilon;

%states
angle = (x(params.ndof)-pi/2);

% Force in the series elastic element
L_m = l0-angle*d;
dL_mdangle = -d;
if strcmp(params.muscletype, 'SEE_linear')
    ksee_lin = ksee*0.005;
    Fsee1 = ksee_lin.*fmax.*(L_m-lce.*lceopt - lsees);
    dFsee1dx = ksee_lin.*fmax.*dL_mdangle;
    dFsee1dlce = diag(-ksee_lin.*fmax.*lceopt);
    
    Fsee = 0.1*fmax.*(L_m-lce.*lceopt-lsees) +  1/2*(Fsee1+sqrt(Fsee1.^2+epsilon^2));
    dFseedFsee1 = 1/2*(1+Fsee1./sqrt(Fsee1.^2+epsilon^2));
    dFseedangle =  0.1*fmax.*dL_mdangle + bsxfun(@times,dFseedFsee1,dFsee1dx);
    dFseedlce = diag(-0.1*fmax.*lceopt) + bsxfun(@times,dFseedFsee1,dFsee1dlce);
    
    stiffness = ksee_lin.*fmax;
elseif strcmp(params.muscletype, 'SEE_quad')
    Fsee = 0.1*fmax.*(L_m-lce.*lceopt-lsees)+(L_m-lce.*lceopt-lsees > 0).*ksee.*fmax.*(L_m-lce.*lceopt-lsees).^2;
    dFseedangle = 0.1*fmax.*dL_mdangle+(L_m-lce.*lceopt-lsees > 0)*2.*ksee.*fmax.*(L_m-lce.*lceopt-lsees).*dL_mdangle;
    dFseedlce = diag(-0.1*fmax.*lceopt-(L_m-lce.*lceopt-lsees>0)*2.*ksee.*fmax.*(L_m-lce.*lceopt-lsees).*lceopt);
    
    stiffness = (L_m-lce.*lceopt-lsees > 0).*ksee.*fmax.*(L_m-lce.*lceopt-lsees)*2;
end
dFseedx = [dFseedangle [0; 0]];


