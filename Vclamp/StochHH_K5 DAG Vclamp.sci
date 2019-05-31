// Potassium channel from original HH model
// Voltage clamp simulations with non-stationary noise analysis
// UNcoupled activation particles (2-state independent particles), 
// Goldwyn et al. (Phys Rev E 83:041908 (2011)) implementation of the
// Diffusion Approximation. Coupled activation particles WITHOUT 
// steady state approximation of variables in the stochastic terms
// See "StochHH_K2 F1 Vclamp noise.sci" for more comments

stacksize('max');
nsim=200;
Tstop=6; dt=0.001; //ms
points = round(Tstop/dt)
NK=300;
//Iamp=10;

Vhold=-90;
Vtest=70;
rand('normal');

tic();

p=1;
Norec = zeros(points,nsim);
//    mrec = zeros(points,nsim);
//    nrec=mrec;
v = Vhold*ones(1,nsim);
an=0.01*(v+55)./(1-exp(-(v+55)/10));
bn=0.125*exp(-(v+65)/80);
N=an./bn;
Kstatesum=(1+N)^4;
n=[ones(1,nsim);4*N;6*N.^2;4*N.^3;N.^4]./(ones(5,1)*Kstatesum);

v = Vtest;
an=0.01*(v+55)./(1-exp(-(v+55)/10));
bn=0.125*exp(-(v+65)/80);

tint=1;
for tt=dt:tint:Tstop
    for t = tt:dt:tt+tint-dt
        Norec(p,:) = n(5,:)*NK;
        p=p+1;

        trans_n=[-(bn+3*an).*n(2,:)+4*an.*n(1,:)+2*bn.*n(3,:);
        -(2*bn+2*an).*n(3,:)+3*an.*n(2,:)+3*bn.*n(4,:);
        -(3*bn+an).*n(4,:)+2*an.*n(3,:)+4*bn.*n(5,:);
        -4*bn.*n(5,:)+an.*n(4,:)];
        
        trans_n=[-(bn+3*an).*n(2,:)+4*an.*n(1,:)+2*bn.*n(3,:);
            -(2*bn+2*an).*n(3,:)+3*an.*n(2,:)+3*bn.*n(4,:);
            -(3*bn+an).*n(4,:)+2*an.*n(3,:)+4*bn.*n(5,:);
            -4*bn.*n(5,:)+an.*n(4,:)];
        
        na=abs(n);
        
        for ii = 1:nsim  //because this algorithm needs a matrix operation (line 90) it has to be done iterating over independent simulations
            Dmtx = [4*an*na(1,ii)+(3*an+bn)*na(2,ii)+2*bn*na(3,ii),-(3*an*na(2,ii)+2*bn*na(3,ii)),0,0;
                -(3*an*na(2,ii)+2*bn*na(3,ii)),3*an*na(2,ii)+2*(an+bn)*na(3,ii)+3*bn*na(4,ii),-(2*an*na(3,ii)+3*bn*na(4,ii)),0;
                0,-(2*an*na(3,ii)+3*bn*na(4,ii)),2*an*na(3,ii)+(an+3*bn)*na(4,ii)+4*bn*na(5,ii),-(an*na(4,ii)+4*bn*na(5,ii));
                0,0,-(an*na(4,ii)+4*bn*na(5,ii)),an*na(4,ii)+4*bn*na(5,ii)]/NK;
            
            Rvec_n(:,ii)=sqrtm(dt*Dmtx)*rand(4,1);
        end
        
        n(2:5,:)=n(2:5,:)+dt*trans_n+Rvec_n;
        n(1,:)=ones(1,nsim)-sum(n(2:5,:),1);

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


