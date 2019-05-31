// Rubinstein model of Ranvier Node's action potentials as in
// Rubinstein JT (1995) Biophys J 68: 779-785.
// Mino H, Rubinstein JT, White JA (2002) Ann Biomed Eng 30: 578-587.
// Bruce IC (2007) Ann Biomed Eng 35: 315-318; 

// Voltage is shifted so that resting voltage = 0 mV (the reversal of the leak)

// This script simulates 1000 sweeps and calculates
// firing efficiency, mean firing time and firing time variance
// 100 (nsim) simultaneous sweeps, 10 times
// Lines 49, 62, 63, 86 and 87 are intended to make a voltage trace
// It is advisable to reduce nsim before uncommenting them
//
// UNcoupled (two-state) activation particles, 3 m and 1 h particle per Na channel
// Difussion Approximation algorithm (F algorithm) with mean values of the variables in the random terms

nsim=100;
gNa=0.02569; //mS/cm2
ENa=144; //mV
Cm=0.0000714; Rm=1953.49;
NNa=1000;
Tstop=1; dt=0.0001; //ms
points = round(Tstop/dt)
Idel=0; Idur=0.1; //ms, ms
threshold=80;

rand('normal');
Eff=[];
meanFT=[];
varFT=[];

tic()

currents=[5:0.1:6.5];
NNa3=3*NNa;
Eff=[];
meanFT=[];
varFT=[];
NaNs=[];
tic()
for curr=1:size(currents,'*')
    Iamp=ones(1,nsim)*currents(curr)/265;
    firetimes=[];
    NNaNs=0;
    for nn=1:10

        p=1;
        //If you want to uncomment the following reduce nsim
        //vrec = zeros(points,nsim);
        v = 0*ones(1,nsim);

        alpha_m=1.872*(v-25.41)./(1-exp((25.41-v)/6.06));
        beta_m=3.973*(21-v)./(1-exp((v-21)/9.41));
        alpha_h=-0.549*(27.74+v)./(1-exp((v+27.74)/9.06));
        beta_h=(22.57)./(1+exp((56-v)/12.5));
        m=ones(1,nsim)./(1+beta_m./alpha_m);
        h=ones(1,nsim)./(1+beta_h./alpha_h);
        firetime=zeros(1,nsim);
        firing=zeros(1,nsim);

        for t = dt:dt:Tstop
            //vrec(p,:) = v;  //uncomment with a reduced nsim
            //p=p+1;
            if or(~firing&v>=threshold) then
                ind=find(v>=threshold&~firing);
                firetime(ind)=t;
                firing(ind)=1;
            end

            alpha_m=1.872*(v-25.41)./(1-exp((25.41-v)/6.06));
            beta_m=3.973*(21-v)./(1-exp((v-21)/9.41));
            alpha_h=-0.549*(27.74+v)./(1-exp((v+27.74)/9.06));
            beta_h=(22.57)./(1+exp((56-v)/12.5));
            SDm = sqrt(2*alpha_m.*beta_m./(alpha_m+beta_m));
            SDh = sqrt(2*alpha_h.*beta_h./(alpha_h+beta_h));    

            Iapp=Iamp*(t>Idel&t<(Idel+Idur));
            Imemb=-Iapp+gNa*m^(3).*h.*(v-ENa)+v./Rm;

            m=m+dt*(alpha_m.*(1-m)-beta_m.*m+rand(1,nsim).*SDm/sqrt(3*NNa*dt));
            h=h+dt*(alpha_h.*(1-h)-beta_h.*h+rand(1,nsim).*SDh/sqrt(NNa*dt));

            v=v-dt*Imemb/Cm;
        end
        printf("round %d Iamp %g\n",nn,currents(curr))
        //clf //uncomment with a reduced nsim
        //plot(dt:dt:Tstop,vrec)
        firetimes=[firetimes firetime];
        NNaNs=NNaNs + sum(1*(isinf(v)|isnan(v)));
    end
    Eff=[Eff;size(find(firetimes<>0),'*')];
    meanFT=[meanFT;mean(firetimes(find(firetimes<>0)))];
    if size(find(firetimes<>0),'*')>1 then
        varFT=[varFT;variance(firetimes(find(firetimes<>0)))];
    else
        varFT=[varFT;0];
    end
    NaNs=[NaNs;NNaNs];
end
realt=toc();
printf("time = %g",realt)

fprintfMat('EffN-Rb2 Fss-N'+string(NNa)+'-dt'+string(dt)+'-'+string(realt)+'s.txt',[currents' Eff meanFT varFT NaNs],'%g\t')

scf(0);
clf
subplot(3,1,1)
plot(currents,Eff)
subplot(3,1,2)
plot(currents,meanFT)
subplot(3,1,3)
plot(currents,varFT)
