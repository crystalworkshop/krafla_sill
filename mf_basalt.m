function CF=mf_basalt(T)
x=T/1000;
 a =       517.9;  
 b =       -1619; 
 c =        1699;
 d =      -597.4;
 CF=1./(1.+exp(a+b*x+c*x.^2+d*x.^3));
end

