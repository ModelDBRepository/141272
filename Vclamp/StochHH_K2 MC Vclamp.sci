// Potassium channel from original HH model
// Voltage clamp simulations with non-stationary noise analysis
// UNcoupled activation particles (2-state independent particles), Diffusion approximation algorithm
// See "StochHH_K5 CWopt Vclamp noise.sci" for more comments

stacksize('max');
nsim=200;  //number of sweeps to be simulated
Tstop=6; dt=0.001; //Total time and dt in ms
points = round(Tstop/dt) //number of points per sweep
NK=300; //number of potassium channels

Vhold=-90;  //voltage at t=0
Vtest=70;
rand('uniform');

tic();

p=1;
Norec = zeros(points,nsim);
v = Vhold*ones(1,nsim);
alpha_n=0.01*(v+55)./(1-exp(-(v+55)/10));
beta_n=0.125*exp(-(v+65)/80);
n=ones(1,nsim)./(1+beta_n./alpha_n);
Nn=round(n*NK*4); //number of ACTIVE 'n' particles (the total is 4*NK)

v = Vtest*ones(1,nsim);
alpha_n=0.01*(v+55)./(1-exp(-(v+55)/10));
beta_n=0.125*exp(-(v+65)/80);
next_evn=-log(1-rand())./(Nn.*beta_n + (NK*4-Nn).*alpha_n);

for t = dt:dt:Tstop
    Norec(p,:)=round(NK*(Nn/(4*NK)).^4);
    p=p+1;
    while or(t>=next_evn)
        //'ii' contains the indices of simulations (sweeps) where a transition has to occur
        ii=find(t>=next_evn); 
        //If rand() < active*beta / (active*beta + inactive*alpha)
        //Then a particle will inactivate; otherwise a particle will activate
        Nn(ii)=Nn(ii)+1-2*(rand(1,size(ii,'*'))<=(Nn(ii).*beta_n(ii)./(Nn(ii).*beta_n(ii)+(NK*4-Nn(ii)).*alpha_n(ii))));
        prev_ev=next_evn(ii);
        next_evn(ii)=prev_ev-log(1-rand())./(Nn(ii).*beta_n(ii) + (NK*4-Nn(ii)).*alpha_n(ii));
    end
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

