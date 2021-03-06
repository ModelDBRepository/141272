% Rubinstein model of Ranvier Node's action potentials as in
% Rubinstein JT (1995) Biophys J 68: 779-785.
% Mino H, Rubinstein JT, White JA (2002) Ann Biomed Eng 30: 578-587.
% Bruce IC (2007) Ann Biomed Eng 35: 315-318; 

% Voltage is shifted so that resting voltage = 0 mV (the reversal of the leak)

% This script simulates 1000 sweeps per stimulus and calculates
% firing efficiency, mean firing time and firing time variance
% 100 (nsim) simultaneous sweeps, 10 times
% Lines 44, 60, 61, 112 and 113 are intended to make a voltage trace
% It is advisable to reduce nsim before uncommenting them
%
% Coupled activation particles, 8-state Na channel
% Diffusion Approximation with steady state values of variables in the random term

nsim=100;
gNa=0.02569*265; %mS/cm2
ENa=144; %mV
Cm=0.0000714*265; Rm=1953.49/265;
NNa=1000;
Tstop=1; dt=0.001; %ms
points = round(Tstop/dt);
Idel=0; Idur=0.1; %ms, ms
threshold=80;
NaNs=0;

currents=5:0.1:6.5;
Eff=[];
meanFT=[];
varFT=[];
NaNs=[];

tic()
for curr=1:length(currents)
    Iamp=ones(1,nsim)*currents(curr);
    firetimes=[];
    NNaNs=0;
    for nn=1:10

        p=1;
        % Uncomment this only with a small nsim
        % vrec = zeros(points,nsim);

        v = 0*ones(1,nsim);
        am=1.872*(v-25.41)./(1-exp((25.41-v)/6.06));
        bm=3.973*(21-v)./(1-exp((v-21)/9.41));
        ah=-0.549*(27.74+v)./(1-exp((v+27.74)/9.06));
        bh=(22.57)./(1+exp((56-v)/12.5));
        M=am./bm;
        H=ah./bh;
        statesum=(1+H).*(1+M).^3;
        s=[ones(1,nsim);3*M;3*M.^2;M.^3;H;3*M.*H;3*M.^(2).*H;M.^(3).*H]./(ones(8,1)*statesum);      

        firetime=zeros(1,nsim);
        firing=zeros(1,nsim);

        for t = dt:dt:Tstop
            %vrec(p,:) = v;  %uncomment with small nsim
            %p=p+1;
            if any(~firing&v>=threshold)
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
            
            % deterministic transition matrix
            transition=[3*am.*s(1,:)-(2*am+bm+ah).*s(2,:)+2*bm.*s(3,:)+bh.*s(6,:);
            2*am.*s(2,:)-(2*bm+am+ah).*s(3,:)+3*bm.*s(4,:)+bh.*s(7,:);
            am.*s(3,:)-(3*bm+ah).*s(4,:)+bh.*s(8,:);
            -(3*am+bh).*s(5,:)+bm.*s(6,:)+ah.*s(1,:);
            3*am.*s(5,:)-(2*am+bm+bh).*s(6,:)+2*bm.*s(7,:)+ah.*s(2,:);
            2*am.*s(6,:)-(2*bm+am+bh).*s(7,:)+3*bm.*s(8,:)+ah.*s(3,:);
            am.*s(7,:)-(3*bm+bh).*s(8,:)+ah.*s(4,:)];

            %steady-state values of the variables
            s_ss=[bm.^(3).*bh;                
            3*am.*bm.^(2).*bh;
            3*am.^(2).*bm.*bh;
            am.^(3).*bh;
            bm.^(3).*ah;
            3*am.*bm.^(2).*ah;
            3*am.^(2).*bm.*ah;
            am.^(3).*ah]./(ones(8,1)*((ah+bh).*(am+bm).^3));

            for ii = 1:nsim  %because this algorithm needs a matrix operation (line 104) it has to be done iterating over independent simulations
                ami=am(ii);bmi=bm(ii);ahi=ah(ii);bhi=bh(ii);
                Dmtx=[3*ami.*s_ss(1,ii)+(2*ami+bmi+ahi).*s_ss(2,ii)+2*bmi.*s_ss(3,ii)+bhi.*s_ss(6,ii),-2*(ami.*s_ss(2,ii)+bmi.*s_ss(3,ii)),0,0,-(ahi.*s_ss(2,ii)+bhi.*s_ss(6,ii)),0,0;
                -2*(ami.*s_ss(2,ii)+bmi.*s_ss(3,ii)),2*ami.*s_ss(2,ii)+(ami+2*bmi+ahi).*s_ss(3,ii)+3*bmi.*s_ss(4,ii)+bhi.*s_ss(7,ii),-(ami.*s_ss(3,ii)+3*bmi.*s_ss(4,ii)),0,0,-(ahi.*s_ss(3,ii)+bhi.*s_ss(7,ii)),0;
                0,-(ami.*s_ss(3,ii)+3*bmi.*s_ss(4,ii)),ami.*s_ss(3,ii)+(3*bmi+ahi).*s_ss(4,ii)+bhi.*s_ss(8,ii),0,0,0,-(ahi*s_ss(4,ii)+bhi.*s_ss(8,ii));
                0,0,0,ahi.*s_ss(1,ii)+(3*ami+bhi).*s_ss(5,ii)+bmi.*s_ss(6,ii),-(3*ami.*s_ss(5,ii)+bmi.*s_ss(6,ii)),0,0;
                -(ahi.*s_ss(2,ii)+bhi.*s_ss(6,ii)),0,0,-(3*ami.*s_ss(5,ii)+bmi.*s_ss(6,ii)),ahi.*s_ss(2,ii)+3*ami.*s_ss(5,ii)+(2*ami+bmi+bhi).*s_ss(6,ii)+2*bmi.*s_ss(7,ii),-2*(ami.*s_ss(6,ii)+bmi.*s_ss(7,ii)),0;
                0,-(ahi.*s_ss(3,ii)+bhi.*s_ss(7,ii)),0,0,-2*(ami.*s_ss(6,ii)+bmi.*s_ss(7,ii)),ahi.*s_ss(3,ii)+2*ami.*s_ss(6,ii)+(ami+2*bmi+bhi).*s_ss(7,ii)+3*bmi.*s_ss(8,ii),-(ami.*s_ss(7,ii)+3*bmi.*s_ss(8,ii));
                0,0,-(ahi*s_ss(4,ii)+bhi.*s_ss(8,ii)),0,0,-(ami.*s_ss(7,ii)+3*bmi.*s_ss(8,ii)),ahi.*s_ss(4,ii)+ami.*s_ss(7,ii)+(3*bmi+bhi).*s_ss(8,ii)]/NNa;

                Rvec(:,ii)=sqrtm(dt*Dmtx)*randn(7,1);
            end

            s(2:8,:)=s(2:8,:)+dt*transition+Rvec;
            s(1,:)=ones(1,nsim)-sum(s(2:8,:),1);

            v=v-dt*Imemb/Cm;
        end
        fprintf('round %d Iamp %g\n',nn,currents(curr))
        % clf %uncomment with small nsim
        % plot(dt:dt:Tstop,vrec)
        firetimes=[firetimes firetime];
        NNaNs=NNaNs + sum(1*(isinf(v)|isnan(v)));
    end
    Eff=[Eff;length(find(firetimes~=0))];
    meanFT=[meanFT;mean(firetimes(firetimes~=0))];
    if length(find(firetimes~=0))>1
        varFT=[varFT;var(firetimes(firetimes~=0))];
    else
        varFT=[varFT;0];
    end
    NaNs=[NaNs;NNaNs];
end

realt=toc();
fprintf('time = %g\n',realt)

csvwrite(sprintf('EffN-Rb8DAGss-N%g-dt%g-%gs.csv',NNa,dt,realt),[currents' Eff meanFT varFT NaNs])

figure(1);
clf
subplot(3,1,1)
plot(currents,Eff)
subplot(3,1,2)
plot(currents,meanFT)
subplot(3,1,3)
plot(currents,varFT)
