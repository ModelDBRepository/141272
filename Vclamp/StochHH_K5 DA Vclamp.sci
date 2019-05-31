// Potassium channel from original HH model
// Voltage clamp simulations with non-stationary noise analysis
// Coupled activation particles (5-state channel), Diffusion approximation algorithm

stacksize('max');
nsim=200;  //number of sweeps to be simulated
Tstop=6; dt=0.001; //Total time and dt in ms
points = round(Tstop/dt) //number of points per sweep
NK=300; //number of potassium channels

Vhold=-90; //voltage for t=0
Vtest=70;
rand('normal');

p=1;
Norec = zeros(points,nsim);

v = Vhold*ones(1,nsim);  //calculus of equilibrium state at t=0
an=0.01*(v+55)./(1-exp(-(v+55)/10));
bn=0.125*exp(-(v+65)/80);
N=an./bn;
Kstatesum=(1+N)^4;
n=[ones(1,nsim);4*N;6*N.^2;4*N.^3;N.^4]./(ones(5,1)*Kstatesum);

//Now we change voltage for the rest of the simulation
v = Vtest*ones(1,nsim);
an=0.01*(v+55)./(1-exp(-(v+55)/10));
bn=0.125*exp(-(v+65)/80);

tic();

tint=1; //period for reporting simulation time (see lines 33 and 61)
for tt=dt:tint:Tstop //Nested FORs are only for the purpose of reporting the time (see line 61)
    for t = tt:dt:tt+tint-dt
        Norec(p,:) = n(5,:)*NK;
        p=p+1;

        trans_n=[-4*an.*n(1,:)+bn.*n(2,:);  //Deterministic part of state transitions
        -(bn+3*an).*n(2,:)+4*an.*n(1,:)+2*bn.*n(3,:);
        -(2*bn+2*an).*n(3,:)+3*an.*n(2,:)+3*bn.*n(4,:);
        -(3*bn+an).*n(4,:)+2*an.*n(3,:)+4*bn.*n(5,:);
        -4*bn.*n(5,:)+an.*n(4,:)];

        na=abs(n); //For the stochastic term we use absolute values of states
        // Stochastic terms R1-R4
        R=rand(4,nsim).*sqrt((dt/NK)*[4*an.*na(1,:)+bn.*na(2,:);
        3*an.*na(2,:)+2*bn.*na(3,:);
        2*an.*na(3,:)+3*bn.*na(4,:);
        an.*na(4,:)+4*bn.*na(5,:)]);

        Wtn=[R(1,:);
        -R(1,:)+R(2,:);
        -R(2,:)+R(3,:);
        -R(3,:)+R(4,:);
        -R(4,:);]
        
        n=n+dt*trans_n+Wtn;
        n(1,:) = ones(1,nsim) - sum(n(2:5,:),1)  //adjusting the sum of states equal to 1

    end
    printf("time %g ms\n",t)
end

time=toc()

scf(0);
clf
plot(dt:dt:Tstop,Norec)

scf(1);
clf
plot(dt:dt:Tstop,[mean(Norec,2),variance(Norec,2)])

scf(2);
clf
plot(mean(Norec,2),variance(Norec,2))


printf("time = %g\n",time);


