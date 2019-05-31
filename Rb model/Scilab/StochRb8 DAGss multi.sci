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
// Goldwin et al. implementation (Phys Rev E 83:041908 (2011))
// Difussion Approximation algorithm with mean values of the variables in the random term

nsim=100;
gNa=0.02569; //mS/cm2
//gNa=0
ENa=144; //mV
Cm=0.0000714; Rm=1953.49;
NNa=1000;
Tstop=1; dt=0.0002; //ms
points = round(Tstop/dt)
Idel=0; Idur=0.1; //�A, ms, ms
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
        //		vrec = zeros(points,nsim);
        v = 0*ones(1,nsim);

        am=1.872*(v-25.41)./(1-exp((25.41-v)/6.06));
        bm=3.973*(21-v)./(1-exp((v-21)/9.41));
        ah=-0.549*(27.74+v)./(1-exp((v+27.74)/9.06));
        bh=(22.57)./(1+exp((56-v)/12.5));
        M=am./bm;
        H=ah./bh;
        statesum=(1+H).*(1+M)^3;
        s=[ones(1,nsim);3*M;3*M.^2;M.^3;H;3*M.*H;3*M.^(2).*H;M^(3).*H]./(ones(8,1)*statesum);  
        firetime=zeros(1,nsim);
        firing=zeros(1,nsim);

        //SDm=0;
        //SDh=0;

        for t = dt:dt:Tstop

            if or(~firing&v>=threshold) then
                ind=find(v>=threshold&~firing);
                firetime(ind)=t;
                firing(ind)=1;
            end

            am=1.872*(v-25.41)./(1-exp((25.41-v)/6.06));
            bm=3.973*(21-v)./(1-exp((v-21)/9.41));
            ah=-0.549*(27.74+v)./(1-exp((v+27.74)/9.06));
            bh=(22.57)./(1+exp((56-v)/12.5));

            Iapp=Iamp*(t>Idel&t<(Idel+Idur));
            Imemb=-Iapp+gNa.*s(8,:).*(v-ENa)+v/Rm;

            transition=[3*am.*s(1,:)-(2*am+bm+ah).*s(2,:)+2*bm.*s(3,:)+bh.*s(6,:);
            2*am.*s(2,:)-(2*bm+am+ah).*s(3,:)+3*bm.*s(4,:)+bh.*s(7,:);
            am.*s(3,:)-(3*bm+ah).*s(4,:)+bh.*s(8,:);
            -(3*am+bh).*s(5,:)+bm.*s(6,:)+ah.*s(1,:);
            3*am.*s(5,:)-(2*am+bm+bh).*s(6,:)+2*bm.*s(7,:)+ah.*s(2,:);
            2*am.*s(6,:)-(2*bm+am+bh).*s(7,:)+3*bm.*s(8,:)+ah.*s(3,:);
            am.*s(7,:)-(3*bm+bh).*s(8,:)+ah.*s(4,:)];

            s_ss=[bm.^(3).*bh;
            3*am.*bm.^(2).*bh;
            3*am.^(2).*bm.*bh;
            am.^(3).*bh;
            bm.^(3).*ah;
            3*am.*bm.^(2).*ah;
            3*am.^(2).*bm.*ah;
            am.^(3).*ah]./(ones(8,1)*((ah+bh).*(am+bm).^3));

            for ii = 1:nsim  //because this algorithm needs a matrix operation (line 109) it has to be done iterating over independent simulations
                ami=am(ii);bmi=bm(ii);ahi=ah(ii);bhi=bh(ii);
                Dmtx=[3*ami.*s_ss(1,ii)+(2*ami+bmi+ahi).*s_ss(2,ii)+2*bmi.*s_ss(3,ii)+bhi.*s_ss(6,ii),-2*(ami.*s_ss(2,ii)+bmi.*s_ss(3,ii)),0,0,-(ahi.*s_ss(2,ii)+bhi.*s_ss(6,ii)),0,0;
                -2*(ami.*s_ss(2,ii)+bmi.*s_ss(3,ii)),2*ami.*s_ss(2,ii)+(ami+2*bmi+ahi).*s_ss(3,ii)+3*bmi.*s_ss(4,ii)+bhi.*s_ss(7,ii),-(ami.*s_ss(3,ii)+3*bmi.*s_ss(4,ii)),0,0,-(ahi.*s_ss(3,ii)+bhi.*s_ss(7,ii)),0;
                0,-(ami.*s_ss(3,ii)+3*bmi.*s_ss(4,ii)),ami.*s_ss(3,ii)+(3*bmi+ahi).*s_ss(4,ii)+bhi.*s_ss(8,ii),0,0,0,-(ahi*s_ss(4,ii)+bhi.*s_ss(8,ii));
                0,0,0,ahi.*s_ss(1,ii)+(3*ami+bhi).*s_ss(5,ii)+bmi.*s_ss(6,ii),-(3*ami.*s_ss(5,ii)+bmi.*s_ss(6,ii)),0,0;
                -(ahi.*s_ss(2,ii)+bhi.*s_ss(6,ii)),0,0,-(3*ami.*s_ss(5,ii)+bmi.*s_ss(6,ii)),ahi.*s_ss(2,ii)+3*ami.*s_ss(5,ii)+(2*ami+bmi+bhi).*s_ss(6,ii)+2*bmi.*s_ss(7,ii),-2*(ami.*s_ss(6,ii)+bmi.*s_ss(7,ii)),0;
                0,-(ahi.*s_ss(3,ii)+bhi.*s_ss(7,ii)),0,0,-2*(ami.*s_ss(6,ii)+bmi.*s_ss(7,ii)),ahi.*s_ss(3,ii)+2*ami.*s_ss(6,ii)+(ami+2*bmi+bhi).*s_ss(7,ii)+3*bmi.*s_ss(8,ii),-(ami.*s_ss(7,ii)+3*bmi.*s_ss(8,ii));
                0,0,-(ahi*s_ss(4,ii)+bhi.*s_ss(8,ii)),0,0,-(ami.*s_ss(7,ii)+3*bmi.*s_ss(8,ii)),ahi.*s_ss(4,ii)+ami.*s_ss(7,ii)+(3*bmi+bhi).*s_ss(8,ii)]/NNa;

                Rvec(:,ii)=sqrtm(dt*Dmtx)*rand(7,1);
            end

            s(2:8,:)=s(2:8,:)+dt*transition+Rvec;

            s(1,:)=ones(1,nsim)-sum(s(2:8,:),'r');

            v=v-dt*Imemb/Cm;
        end
        printf("round %d Iamp %g\n",nn,currents(curr))
        //		clf
        //		plot(dt:dt:Tstop,vrec)
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

write_csv([currents' Eff meanFT varFT NaNs],sprintf('EffN-Rb8DAGss-N%g-dt%g-%gs.csv',NNa,dt,realt),',',".")

//scf(0);
//clf
//subplot(3,1,1)
//plot(currents,Eff)
//subplot(3,1,2)
//plot(currents,meanFT)
//subplot(3,1,3)
//plot(currents,varFT)
