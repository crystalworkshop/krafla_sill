function c=heatcap(x,T,par)
if par.bas==0

%dXdt for rhyolite Pivinskii
% t1 = T .^ 2;
% t11 = exp(0.961026371399999562e3 - 0.1866187556e-5 .* t1 .* T + t1 * 0.447948339799999953e-2 + T * (-0.359050896099999894e1));
% t14 = (0.1e1 + t11) .^ 2;
% dbetadtrh = 0.5598562668e-5 ./ t14 .* t11 .* (t1 - 0.1600226224e4 * T + 0.6413269216e6);
t1 = T .^ 2;
t9 = exp(0.3043e4 - 0.471e-5 * t1 .* T + t1 * 0.1224e-1 + T .* (-0.105e2));
t12 = (0.1e1 + t9) .^ 2;
dbetadtrh = 0.1412700000e-4 ./ t12 * t9 .* (t1 - 0.1727882778e4 * T + 0.7469384866e6);
c = (par.Cp_rh+dbetadtrh*par.Lh)*par.rho_rh;

else

%dXdt for rhyolite Piwinsky
% t1 = T .^ 2;
% t11 = exp(0.961026371399999562e3 - 0.1866187556e-5 .* t1 .* T + t1 * 0.447948339799999953e-2 + T * (-0.359050896099999894e1));
% t14 = (0.1e1 + t11) .^ 2;
% dbetadtrh = 0.5598562668e-5 ./ t14 .* t11 .* (t1 - 0.1600226224e4 * T + 0.6413269216e6);
dT=par.dT;
t1 = (T-dT) .^ 2;
t9 = exp(0.3043e4 - 0.4709000000e-5 .* t1 .* (T-dT) + t1 * 0.122048999999999994e-1 + (T-dT) .* (-0.105519999999999996e2));
t12 = (0.1e1 + t9) .^ 2;
dbetadtrh = 0.1412700000e-4 ./ t12 .* t9 .* (t1 - 0.1727882778e4 * (T-dT) + 0.7469384866e6);
  
%dXdt for basalt
t1 = T.^ 2;
t9 = exp(0.517900000000004638e3 - 0.5974000000e-6 * t1 .* T + t1 * 0.169899999999999999e-2 + T * (-0.161900000000000799e1));
t12 = (0.1e1 + t9) .^ 2;
dbetadtbas = 0.1792200000e-5 ./ t12 .* t9 .* (t1 - 0.1895993751e4 * T + 0.9033590000e6);
    
c(x<0)  = (par.Cp_bas+dbetadtbas(x<0)*par.Lh)*par.rho_bas;
c(x>=0) = (par.Cp_rh +dbetadtrh(x>=0)*par.Lh)*par.rho_rh;
% c=par.Cp_rh*par.rho_rh*ones(size(T));
end