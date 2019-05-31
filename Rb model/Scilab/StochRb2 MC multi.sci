// Rubinstein model of Ranvier Node's action potentials as in
// Rubinstein JT (1995) Biophys J 68: 779-785.
// Mino H, Rubinstein JT, White JA (2002) Ann Biomed Eng 30: 578-587.
// Bruce IC (2007) Ann Biomed Eng 35: 315-318; 

// Voltage is shifted so that resting voltage = 0 mV (the reversal of the leak)

// This script simulates 1000 sweeps and calculates
// firing efficiency, mean firing time and firing time variance
// 100 (nsim) simultaneous sweeps, 10 times
// Lines 35, 51, 52, 84 and 85 are intended to make a voltage trace
// It is advisable to reduce nsim before uncommenting them
//
// UNcoupled (two-state) activation particles, 3 m and 1 h particle per Na channel
// Markov Chain modeling (Chow&White algorithm)

nsim=100;
gNa=0.02569*265; //mS/cm2
ENa=144; //mV
Cm=0.0000714*265; Rm=1953.49/265;
NNa=1000; //Number of sodium channels
Tstop=1; dt=0.001; //ms
Idel=0; Idur=0.1; //µA, ms, ms
threshold=80;
points = round(Tstop/dt)

rand('uniform');

currents=[5:0.1:6.5];
NNa3=3*NNa;
Eff=[];
meanFT=[];
varFT=[];
NaNs=[];
tic()
for curr=1:size(currents,'*')
    Iamp=ones(1,nsim)*currents(curr);
    firetimes=[];
    NNaNs=0;
    for nn=1:10
        p=1;
        //If you want to uncomment the following reduce nsim
        //vrec = zeros(points,nsim);
        v = zeros(1,nsim);
        alpha_m=1.872*(v-25.41)./(1-exp((25.41-v)/6.06));
        beta_m=3.973*(21-v)./(1-exp((v-21)/9.41));
        beta_h=(22.57)./(1+exp((56-v)/12.5));
        alpha_h=-0.549*(27.74+v)./(1-exp((v+27.74)/9.06));        
        m=ones(1,nsim)./(1+beta_m./alpha_m);
        Nm=round(m*NNa*3);
        h=ones(1,nsim)./(1+beta_h./alpha_h);
        Nh=round(h*NNa);
        next_evh=-log(rand())./(Nh.*beta_h + (NNa-Nh).*alpha_h); //next 'h' event
        next_evm=-log(rand())./(Nm.*beta_m + (NNa3-Nm).*alpha_m); //next 'm' event
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

            Iapp=Iamp*(t>Idel&t<(Idel+Idur));
            m=Nm/(NNa3); h=Nh/NNa;
            Imemb=-Iapp+gNa*m^(3).*h.*(v-ENa)+v./Rm;

            while or(t>=next_evm)  //activation or deactivation of a single m particle
                ii=find(t>=next_evm); // ii are the simulations in which an m transition is going to occur
                Nm(ii)=Nm(ii)+1-2*(rand(1,size(ii,'*'))<=(Nm(ii).*beta_m(ii)./(Nm(ii).*beta_m(ii)+(NNa3-Nm(ii)).*alpha_m(ii))));
                prev_ev=next_evm(ii);  
                alpha_m=1.872*(v-25.41)./(1-exp((25.41-v)/6.06));
                beta_m=3.973*(21-v)./(1-exp((v-21)/9.41));
                next_evm(ii)=prev_ev-log(rand())./(Nm(ii).*beta_m(ii) + (NNa3-Nm(ii)).*alpha_m(ii));
            end

            while or(t>=next_evh)    //activation or deactivation of a single h particle
                ii=find(t>=next_evh); // ii are the simulations in which an h transition is going to occur
                Nh(ii)=Nh(ii)+1-2*(rand(size(ii,'*'))<=(Nh(ii).*beta_h(ii)./(Nh(ii).*beta_h(ii)+(NNa-Nh(ii)).*alpha_h(ii))));
                prev_ev=next_evh(ii);
                beta_h=(22.57)./(1+exp((56-v)/12.5));
                alpha_h=-0.549*(27.74+v)./(1-exp((v+27.74)/9.06));        
                next_evh(ii)=prev_ev-log(rand())./(Nh(ii).*beta_h(ii) + (NNa-Nh(ii)).*alpha_h(ii));
            end    
            v=v-dt*Imemb/Cm;
        end
        printf("round %d Iamp %g\n",nn,currents(curr))
        //clf  //uncomment with a reduced nsim
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
printf("time = %g\n",realt)
fprintfMat('EffN-Rb2 CW-N'+string(NNa)+'-dt'+string(dt)+'-'+string(realt)+'s.txt',[currents' Eff meanFT varFT NaNs],'%g\t')

scf(0);
clf
subplot(3,1,1)
plot(currents,Eff)
subplot(3,1,2)
plot(currents,meanFT)
subplot(3,1,3)
plot(currents,varFT)
