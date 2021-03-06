// Potassium channel from original HH model
// Voltage clamp simulations with non-stationary noise analysis
// UNcoupled activation particles (2-state independent particles), Diffusion approximation algorithm
// Steady state approximation of variables in the stochastic terms
// See "StochHH_K2 F1 Vclamp noise.sci" for more comments

stacksize('max');
nsim=200;  //number of sweeps to be simulated
Tstop=6; dt=0.001; //Total time and dt in ms
points = round(Tstop/dt) //number of points per sweep
NK=300; //number of potassium channels

Vhold=-90;
Vtest=70;
rand('normal')

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

tic()

for t = dt:dt:Tstop
    Norec(p,:) = NK*n.^4;
    p=p+1;
	SDn = sqrt(abs(2*alpha_n.*beta_n)/((alpha_n+beta_n)*dt*NK*4));
    n=n+dt*(alpha_n.*(1-n)-beta_n.*n+rand(1,nsim).*SDn);
    
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

