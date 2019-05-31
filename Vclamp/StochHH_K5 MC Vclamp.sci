// Potassium channel from original HH model
// Voltage clamp simulations with non-stationary noise analysis
// Coupled activation particles (5-state channel), Markov Chain modeling

stacksize('max')
nsim=200;  //number of sweeps to be simulated
Tstop=6; dt=0.001; //Total time and dt in ms
points = round(Tstop/dt) //number of points per sweep
NK=300; //number of potassium channels

Vhold=-90; //voltage for t=0
Vtest=70;
rand('uniform');

//K_trans will have one column per possible transition (8)
K_trans=[-1 1 0 0 0
1 -1 0 0 0
0 -1 1 0 0
0 1 -1 0 0
0 0 -1 1 0
0 0 1 -1 0
0 0 0 -1 1
0 0 0 1 -1]';
xx=zeros(1,nsim);

p=1;
Norec = zeros(points,nsim);

v = Vhold*ones(1,nsim); //calculus of equilibrium state at t=0
alpha_n=0.01*(v+55)./(1-exp(-(v+55)/10));
beta_n=0.125*exp(-(v+65)/80);
N=alpha_n./beta_n;
Kstatesum=(1+N)^4;
Kstates=round(NK*[ones(1,nsim);4*N;6*N.^2;4*N.^3;N.^4]./(ones(5,1)*Kstatesum));

//Now we change voltage for the rest of the simulation
v = Vtest*ones(1,nsim);
alpha_n=0.01*(v+55)./(1-exp(-(v+55)/10));
beta_n=0.125*exp(-(v+65)/80);

Krates=[4*alpha_n.*Kstates(1,:)
beta_n.*Kstates(2,:)
3*alpha_n.*Kstates(2,:)
2*beta_n.*Kstates(3,:)
2*alpha_n.*Kstates(3,:)
3*beta_n.*Kstates(4,:)
alpha_n.*Kstates(4,:)
4*beta_n.*Kstates(5,:)];

next_evK=-log(rand(1,nsim))./sum(Krates,'r'); //Time for the next transition (one per sweep)
tint = 1; //period for reporting simulation time (see lines 54 and 85)
tic();

for tt=dt:tint:Tstop  //Nested FORs are only for the purpose of reporting the time (see line 85)
    for t = tt:dt:tt+tint-dt

        Norec(p,:) = Kstates(5,:); //this is the number of open channels at time t (position p of Norec)
        p=p+1;

        while or(t>=next_evK) 
            ii=find(t>=next_evK); //ii = in which simulations (sweeps) a transition is going to occur
            dist=cumsum(Krates./(ones(8,1)*sum(Krates,'r')),'r'); //Cummulative probabilities matrix
            ev=rand(1,nsim); 
            for a=ii         
                xx(a)=min(find(ev(a)<dist(:,a))); //choosing a transition (out of 8) for the 'ii' simulation
            end

            Kstates(:,ii)=Kstates(:,ii)+K_trans(:,xx(ii));  

            Krates=[4*alpha_n.*Kstates(1,:)
            beta_n.*Kstates(2,:)
            3*alpha_n.*Kstates(2,:)
            2*beta_n.*Kstates(3,:)
            2*alpha_n.*Kstates(3,:)
            3*beta_n.*Kstates(4,:)
            alpha_n.*Kstates(4,:)
            4*beta_n.*Kstates(5,:)];

            prev_ev=next_evK(ii); //next event is calculated starting from the last
                                  //this and the use of while instead of for, allows more than one transition
                                  // in one time step 
            next_evK(ii)=prev_ev-log(rand(1,size(ii,'*')))./sum(Krates(:,ii),'r');
        end
    end
    printf("time %g ms\n",t);
end

time=toc()
printf("time = %g\n",time);

scf(0);
clf
plot(dt:dt:Tstop,Norec) //comment this if you have low ram

scf(1);
clf
plot(dt:dt:Tstop,[mean(Norec,2),variance(Norec,2)])

scf(2);
clf
plot(mean(Norec,2),variance(Norec,2))



