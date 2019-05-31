// Potassium channel from original HH model
// Voltage clamp simulations with non-stationary noise analysis
// UNcoupled activation particles (2-state independent particles), 
// Linaro et al. (Plos Comput Biol 7:e1001102 (2011)) algorithm for channel noise
// Steady state approximation of variables in the stochastic terms
// See "StochHH_K2 F1 Vclamp noise.sci" for more comments

stacksize('max')
nsim=200;
NK=300;
Tstop=6; dt=0.001; //ms
//Iamp=10;

Vhold=-90;
Vtest=70;
rand('normal')
tic()

points = round(Tstop/dt)
p=1;
Norec = zeros(points,nsim);
v = Vhold*ones(1,nsim);
alpha_n=0.01*(v+55)./(1-exp(-(v+55)/10));
beta_n=0.125*exp(-(v+65)/80);
n=ones(1,nsim)./(1+beta_n./alpha_n);    

v = Vtest*ones(1,nsim);
alpha_n=0.01*(v+55)./(1-exp(-(v+55)/10));
beta_n=0.125*exp(-(v+65)/80);
nss=alpha_n./(alpha_n+beta_n);
tn=(1)./(alpha_n+beta_n);
ChiK=0;

for t = dt:dt:Tstop
    Norec(p,:) = NK*(n.^4 +sum(ChiK,1));
    p=p+1;
       
    n=n+dt*(alpha_n.*(1-n)-beta_n.*n);
    Sign = sqrt([4*nss.^7 .*(1-nss);
                 6*nss.^6 .*(1-nss)^2;
                 4*nss.^5 .*(1-nss)^3;
                 nss.^4 .*(1-nss).^4]./NK);
    Taun = [tn;
            tn/2;
            tn/3;
            tn/4];
            
    ChiK = ChiK + dt*(-ChiK + Sign.*sqrt(2*Taun/dt).*rand(4,nsim))./Taun;
    
end

time=toc()
printf("time = %g\n",time);

scf(0);
clf
plot(dt:dt:Tstop,Norec)

scf(1);
clf
plot(dt:dt:Tstop,[mean(Norec,2),variance(Norec,2)])

scf(2);
clf
plot(mean(Norec,2),variance(Norec,2))

