// Potassium channel from original HH model
// Voltage clamp simulations with non-stationary noise analysis
// UNcoupled activation particles (2-state independent particles), 
// Goldwyn et al. (Phys Rev E 83:041908 (2011)) implementation of the
// Diffusion Approximation. Coupled activation particles with 
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

nss=[bn.^4;
     4*an.*bn.^3;
     6*an.^2 .*bn.^2;
     4*an.^3 .*bn;
     an.^4]./((an + bn).^4);
     
Dmtx = [4*an*nss(1)+(3*an+bn)*nss(2)+2*bn*nss(3),-(3*an*nss(2)+2*bn*nss(3)),0,0;
        -(3*an*nss(2)+2*bn*nss(3)),3*an*nss(2)+2*(an+bn)*nss(3)+3*bn*nss(4),-(2*an*nss(3)+3*bn*nss(4)),0;
        0,-(2*an*nss(3)+3*bn*nss(4)),2*an*nss(3)+(an+3*bn)*nss(4)+4*bn*nss(5),-(an*nss(4)+4*bn*nss(5));
        0,0,-(an*nss(4)+4*bn*nss(5)),an*nss(4)+4*bn*nss(5)]/NK;
        
Smtx=sqrtm(Dmtx*dt);


tint=1;
for tt=dt:tint:Tstop
    for t = tt:dt:tt+tint-dt
        Norec(p,:) = n(5,:)*NK;
        p=p+1;

        trans_n=[-(bn+3*an).*n(2,:)+4*an.*n(1,:)+2*bn.*n(3,:);
        -(2*bn+2*an).*n(3,:)+3*an.*n(2,:)+3*bn.*n(4,:);
        -(3*bn+an).*n(4,:)+2*an.*n(3,:)+4*bn.*n(5,:);
        -4*bn.*n(5,:)+an.*n(4,:)];
        
        n(2:5,:)=n(2:5,:)+dt*trans_n+Smtx*rand(4,nsim);
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


