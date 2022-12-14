function k=therm_cond(x,T,par)

kmin=1;
kmax=100;
% x=xx.*ones(size(T));
mfcrit=0.95;
if par.bas==1

    magrock=find(x>0,1);
%in basalt
    Tx=mean(T(:,1:magrock));
    C=heatcap(x(1),Tx,par);
    kappa=par.lam/C;
    bx=1-mf_basalt(Tx);
    Cx=0.02;
    mu=vism(Cx,Tx+273.15,bx);
    dT=T(:,1)-T(:,magrock);
    Ra=par.rho_bas*9.81*par.beta_bas*dT*(par.width/2)^3/kappa/mu;
    k_bas=0.2*Ra.^(1/3);
%in rhyolite
    molt_rh=find(mf_rh(T)<mfcrit,1);
    if isempty(molt_rh), molt_rh=magrock+1; end
    if molt_rh > magrock+1
        Tx=mean(T(:,magrock+1:molt_rh))+273.15;
    else
        molt_rh = magrock+1;
        Tx=T(:,magrock+1)+273.15;
    end
    C=heatcap(x(end),Tx-273.15,par);
    kappa=par.lam/C;
    bx=1-mf_rh(Tx);
    Cx=0.02;
    mu=vism(Cx,Tx+273.15,bx);
    dT=T(:,magrock+1)-T(:,molt_rh);
    h_rh=x(molt_rh)- x(magrock);
    Ra=par.rho_rh*9.81*par.beta_rh*dT*h_rh^3/kappa/mu;
    k_rh=max(par.lam,0.2*Ra.^(1/3));

    k(:,1:magrock)=k_bas;
    k(:,magrock+1:molt_rh)=k_rh;
    k(molt_rh+1:size(T,2))=par.lam;
    k(2:end)=(k(1:end-1)+k(2:end))/2;
else
    magrock=find(x>0,1);
    Tx=mean(T(:,1:magrock))+273.15;
    kappa=par.lam./heatcap(x,Tx-273.15,par);

    % viscosity
    bx=1-mean(mf_rh(T(:,1:magrock)));
    Cx=0.02;
    mu=vism(Cx,Tx,bx);
    dT=T(:,1)-T(:,magrock+1);
    Ra=par.rho_rh*9.81*par.beta_rh*dT*(par.width/2)^3/kappa/mu;
    kmax=0.2*Ra.^(1/3);
    %     k = (kmax-kmin)*(2*atan((500 * ())) / pi+1)/2 + kmin;
    k = (kmax-kmin)*(2*atan((1000 * (mf_rock(T)-mfcrit))) / pi+1)/2 + kmin;
end

% k(:,x>=0) = (kmax-kmin) *(2*atan((500 * (mf_rock(T(:,x>=0))-mfcrit))) / pi+1)/2 + kmin;
%     k(:,x<0)  = (kmax-kmin) *(2*atan((500 * (mf_basalt(T(:,x<0))-mfcrit))) / pi+1)/2 + kmin;

% k(x>=0)=kmin;
%
% k(x<0) = -0.100e3 * atan((100 * (x(x<0)))) / pi + 0.51e2;
% k=ones(size(x));
end
function vm=vism(cx,Tx,bx) %Hess and Dingwell
t1=cx*100;
vmelt=10.^(-0.3545e1+0.833.*t1+(9601-2368.*t1)./(Tx-195.7-32.25.*t1));
aa1= 0.9999e0;
aa2= 1.53e0;
aa3= 1.91e0;
aa4=13.0e0-aa3;
Fbx=erf(0.88622796e0/aa1.*(aa2*bx)*(1+(aa2*bx).^aa3));
teta2 = (1+(aa2*bx).^aa4).*(1-aa1*Fbx).^(-2.5/aa2);
vm=vmelt.*teta2;
end
