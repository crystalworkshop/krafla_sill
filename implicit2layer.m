clear all
close all
nx = 500;

par.width=150; %widt of the sill
par.Cp_bas=1200; % heat capacity of basalt
par.Cp_rh =1200; % heat capacity of rhyolite
par.lam     =1;  % thermal conductibity 
par.bas     =1; % 1- basalt 0 -rhyolite
par.rho_rh  =2300; %density of basalt
par.rho_bas =2800; %density of rhyolite
par.beta_bas=5e-5; %thermal compressibility
par.beta_rh =5e-5; %thermal compressibility
par.Lh      =3.5e5;%Latent heat of crystallization
par.Tinit=1100;    %initial temperatire
par.dT=50;

x=linspace(-par.width/2,par.width/2+100,5000); %mesh 
nt= 500;
tyear=3600*24*365;
a=-tyear*35; %final time
b=0;
tt=1:nt;
t=a+2*(b-a)*(tt-1)/(nt-1)-(b-a)*((tt-1)/(nt-1)).^2;
t=-t(nt:-1:1); % time vector
dt=diff(t);
n=numel(x);
dx(2:n)=diff(x);
dx(1)=dx(2);

dxi=(dx(1:n-1)+dx(2:n))/2;
dxi(n)=dxi(n-1);
T0= icfun ( x,par );
A=zeros(1,n);B=zeros(1,n);C=zeros(1,n);F=zeros(1,n);
u=zeros(nt,n);
u(1,:)=T0;
for i=2:nt
    eps=1e10;
    Ti=T0;
    it=1;
    while eps>1e-6

        c=heatcap(x,Ti,par);
        k=therm_cond(x,Ti,par);

        B(1)=1;
        C(1)=1;
        F(1)=0;

        A(2:n-1)=dt(i-1)*(1./k(1:n-2)+1./k(2:n-1)).^(-1)./dxi(1:n-2);
        B(2:n-1)=dt(i-1)*(1./k(2:n-1)+1./k(3:n)).^(-1)./dxi(2:n-1);
        C(2:n-1)=A(2:n-1)+B(2:n-1)+c(2:n-1).*dx(2:n-1)/2;
        F(2:n-1)=c(2:n-1).*dx(2:n-1).*T0(2:n-1)/2;

        A(n)=0;
        C(n)=1;
        F(n)=T0(n);

        T=trid(A,B,C,F,n);
        eps=sum(((T-Ti)/500).^2)/n;
        Ti=T;
        it=it+1;
        if it>1000, break; end
    end
    disp(eps)
    T0=T;
    u(i,:)=T;
    if mod(i,10)==0, plot(x,T,'-',LineWidth=1); end
    hold on
    pause(1e-10)
end

%% Plot results
par.tyear=tyear;
figure ( 2 )
surf ( x, t/tyear, u, 'EdgeColor', 'None' );
xlabel ( 'distance, m' )
ylabel ( ' time, years' );
zlabel ( 'Temperature, ^oC' );

figure ( 3 )
k=u;
Cp=u;
xs=t;
for i=1:size(u,1)-1
    nexttile(1)
    Tx=u(i,:);
    k(i,:)=therm_cond(x,Tx,par);
    Cp(i,:)=heatcap(x,Tx,par);
    if mod(i,10) ==0, plot(x,k(i,:)); end
    xlim([-10 50])
    xlabel ( 'distance from init. contact' )
    ylabel ( 'Thermal conductivity' );
    hold on
    nexttile(2)
    if mod(i,10) ==0, plot(x,Cp(i,:)); end
    xlim([-10 50])
    xlabel ( 'distance from init. contact' )
    ylabel ( 'Heat capacity' );
    hold on
   
    ik=find(k(i,:)>5,1,"last");
    if isempty(ik)
        xs(i)=nan;
    else
        xs(i)=x(ik);
    end
    hold on
end
xs(end)=xs(end-1);
nexttile(3)
plot(t/tyear,xs,'-',LineWidth=2)
xlabel ( 'time, year' )
ylabel ( 'melting front position' );


figure(4)
nexttile(1)
imag=find(x>0,1);
ibas=find(x<=-10,1,"last");
irhy=find(x>=30,1);

indbas=[1,imag:-10:ibas];
indrhy=imag+1:10:irhy;
Tbas=u(:,indbas);
Trh =u(:,indrhy );
plot(t/tyear,Tbas,t/tyear,Trh,'--')
ylim([800 par.Tinit])
ind=[indbas,indrhy];
T=[Tbas,Trh];
rn=max(1,fix(rand(size(ind))*nt));
for i=1:numel(ind)
    if T(rn(i),i) > 800, text(t(rn(i))/tyear,T(rn(i),i),num2str(x(ind(i)), '%4.1f'  ),"FontSize",14,"FontWeight","bold"); end
end
xlabel('time,years')
ylabel('Temperature,^oC')
nexttile(2)
basrock=find(x >= par.width/2,1);
if par.bas==1
    mfbas=mf_basalt(Tbas);
else
    mfbas=mf_rh(Tbas);
end
mfrh=mf_rh(Trh);

plot(t/tyear,mfbas, t/tyear,mfrh,'--')
mf=[mfbas,mfrh];
for i=1:numel(ind)
   if mf(rn(i),i) > 0.05,  text(t(rn(i))/tyear,mf(rn(i),i),num2str(x(ind(i)), '%4.1f'),"FontSize",14,"FontWeight","bold"); end
end
xlabel('time,years')
ylabel('Melt fraction')

T=u;
save Thist.mat t T x ind


function u0 = icfun ( x,par )
u0 (x>0) = 400;
u0( x<=0 ) = par.Tinit;
end

function V=trid(A,B,C,F,nx)
persistent alpha1 beta1

if isempty(alpha1)
    alpha1=zeros(1,nx);
    beta1=zeros(1,nx);
end

%Forward Thomas path
alpha1(1)=B(1)/C(1);
beta1(1)= F(1)/C(1);
for i=1:nx-1
    zn=1./(C(i)-A(i)*alpha1(i));
    alpha1(i+1)=B(i)*zn;
    beta1(i+1)=(F(i)+A(i)*beta1(i))*zn;
end
%Backward Thomas path
V(nx)=(F(nx)+A(nx)*beta1(nx))/(C(nx)-A(nx)*alpha1(nx));
for i=nx-1:-1:1
    V(i)=V(i+1)*alpha1(i+1)+beta1(i+1);
end

end
