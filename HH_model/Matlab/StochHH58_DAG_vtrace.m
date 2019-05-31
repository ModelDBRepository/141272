% Stochastic Hodgkin and Huxley model
% Voltage is shifted from original model to agree with current conventions (Vext = 0)
%
% This script will simulate 10 parallel simulations of 50 ms and output the
% voltage traces.
%
% Diffusion approximation algorithm, as described by Fox 1997 (Biophys J.)
% and implemented by Goldywn 2011 (Phys Rev E). The square root of  the
% diffusion matrix is calculated at each time step WITHOUT the steady state
% approximation.

nsim=10;
gK=36; gNa=120; gL=0.3; %mS/cm2
EK=-77; EL=-54.4; ENa=50; %mV
Tstop=50; dt=0.005; %ms
points = round(Tstop/dt);
NNa=6000;
NK=NNa*0.3;%NK=900;

threshold=-10;
Eff=[];

tic();

firetimes=[];
xx=zeros(1,nsim);

vrec = zeros(points,nsim);
p=1;
v = -65*ones(1,nsim);
an=0.01*(v+55)./(1-exp(-(v+55)/10));
bn=0.125*exp(-(v+65)/80);
am=0.1*(v+40)./(1-exp(-(v+40)/10));
bm=4*exp(-(v+65)/18);
bh=(1)./(1+exp(-(v+35)/10));
ah=0.07*exp(-(v+65)/20);
M=am./bm;
H=ah./bh;
Nastatesum=(1+H).*(1+M).^3;
m=[ones(1,nsim);3*M;3*M.^2;M.^3;H;3*M.*H;3*M.^(2).*H;M.^(3).*H]./(ones(8,1)*Nastatesum);

N=an./bn;
Kstatesum=(1+N).^4;
n=[ones(1,nsim);4*N;6*N.^2;4*N.^3;N.^4]./(ones(5,1)*Kstatesum);

spikes=[];
firing=zeros(1,nsim);
tint=5;
for tt=dt:tint:Tstop
    for t = tt:dt:tt+tint-dt
        vrec(p,:) = v;
        p=p+1;
        if any(v>threshold&~firing)
            ind=find(v>threshold&~firing);
            for a=ind;spikes=[spikes;[a t]];end;
            firing(ind)=1;
        end
        
        if any(v<=threshold & firing)
            firing(v<=threshold) = 0;
        end
        
        Imemb=gK*n(5,:).*(v-EK)+gNa*m(8,:).*(v-ENa)+gL*(v-EL);
        
        an=0.01*(v+55)./(1-exp(-(v+55)/10));
        bn=0.125*exp(-(v+65)/80);
        am=0.1*(v+40)./(1-exp(-(v+40)/10));
        bm=4*exp(-(v+65)/18);
        bh=(1)./(1+exp(-(v+35)/10));
        ah=0.07*exp(-(v+65)/20);
        
        trans_m=[3*am.*m(1,:)-(2*am+bm+ah).*m(2,:)+2*bm.*m(3,:)+bh.*m(6,:);
            2*am.*m(2,:)-(2*bm+am+ah).*m(3,:)+3*bm.*m(4,:)+bh.*m(7,:);
            am.*m(3,:)-(3*bm+ah).*m(4,:)+bh.*m(8,:);
            -(3*am+bh).*m(5,:)+bm.*m(6,:)+ah.*m(1,:);
            3*am.*m(5,:)-(2*am+bm+bh).*m(6,:)+2*bm.*m(7,:)+ah.*m(2,:);
            2*am.*m(6,:)-(2*bm+am+bh).*m(7,:)+3*bm.*m(8,:)+ah.*m(3,:);
            am.*m(7,:)-(3*bm+bh).*m(8,:)+ah.*m(4,:)];
        
        ma=abs(m);
        
        for ii = 1:nsim  %because this algorithm needs a matrix operation (line 90) it has to be done iterating over independent simulations
            ami=am(ii);bmi=bm(ii);ahi=ah(ii);bhi=bh(ii);
            Dmtx=[3*ami.*ma(1,ii)+(2*ami+bmi+ahi).*ma(2,ii)+2*bmi.*ma(3,ii)+bhi.*ma(6,ii),-2*(ami.*ma(2,ii)+bmi.*ma(3,ii)),0,0,-(ahi.*ma(2,ii)+bhi.*ma(6,ii)),0,0;
                -2*(ami.*ma(2,ii)+bmi.*ma(3,ii)),2*ami.*ma(2,ii)+(ami+2*bmi+ahi).*ma(3,ii)+3*bmi.*ma(4,ii)+bhi.*ma(7,ii),-(ami.*ma(3,ii)+3*bmi.*ma(4,ii)),0,0,-(ahi.*ma(3,ii)+bhi.*ma(7,ii)),0;
                0,-(ami.*ma(3,ii)+3*bmi.*ma(4,ii)),ami.*ma(3,ii)+(3*bmi+ahi).*ma(4,ii)+bhi.*ma(8,ii),0,0,0,-(ahi*ma(4,ii)+bhi.*ma(8,ii));
                0,0,0,ahi.*ma(1,ii)+(3*ami+bhi).*ma(5,ii)+bmi.*ma(6,ii),-(3*ami.*ma(5,ii)+bmi.*ma(6,ii)),0,0;
                -(ahi.*ma(2,ii)+bhi.*ma(6,ii)),0,0,-(3*ami.*ma(5,ii)+bmi.*ma(6,ii)),ahi.*ma(2,ii)+3*ami.*ma(5,ii)+(2*ami+bmi+bhi).*ma(6,ii)+2*bmi.*ma(7,ii),-2*(ami.*ma(6,ii)+bmi.*ma(7,ii)),0;
                0,-(ahi.*ma(3,ii)+bhi.*ma(7,ii)),0,0,-2*(ami.*ma(6,ii)+bmi.*ma(7,ii)),ahi.*ma(3,ii)+2*ami.*ma(6,ii)+(ami+2*bmi+bhi).*ma(7,ii)+3*bmi.*ma(8,ii),-(ami.*ma(7,ii)+3*bmi.*ma(8,ii));
                0,0,-(ahi*ma(4,ii)+bhi.*ma(8,ii)),0,0,-(ami.*ma(7,ii)+3*bmi.*ma(8,ii)),ahi.*ma(4,ii)+ami.*ma(7,ii)+(3*bmi+bhi).*ma(8,ii)]/NNa;
            
            Rvec(:,ii)=sqrtm(dt*Dmtx)*randn(7,1);
        end
        
        m(2:8,:)=m(2:8,:)+dt*trans_m+Rvec;
        m(1,:)=ones(1,nsim)-sum(m(2:8,:),1);
        
        
        trans_n=[-(bn+3*an).*n(2,:)+4*an.*n(1,:)+2*bn.*n(3,:);
            -(2*bn+2*an).*n(3,:)+3*an.*n(2,:)+3*bn.*n(4,:);
            -(3*bn+an).*n(4,:)+2*an.*n(3,:)+4*bn.*n(5,:);
            -4*bn.*n(5,:)+an.*n(4,:)];
        
        na=abs(n);
        
        for ii = 1:nsim  %because this algorithm needs a matrix operation (line 90) it has to be done iterating over independent simulations
            ani=an(ii);bni=bn(ii);
            Dmtx = [4*ani*na(1,ii)+(3*ani+bni)*na(2,ii)+2*bni*na(3,ii),-(3*ani*na(2,ii)+2*bni*na(3,ii)),0,0;
                -(3*ani*na(2,ii)+2*bni*na(3,ii)),3*ani*na(2,ii)+2*(ani+bni)*na(3,ii)+3*bni*na(4,ii),-(2*ani*na(3,ii)+3*bni*na(4,ii)),0;
                0,-(2*ani*na(3,ii)+3*bni*na(4,ii)),2*ani*na(3,ii)+(ani+3*bni)*na(4,ii)+4*bni*na(5,ii),-(ani*na(4,ii)+4*bni*na(5,ii));
                0,0,-(ani*na(4,ii)+4*bni*na(5,ii)),ani*na(4,ii)+4*bni*na(5,ii)]/NK;
            
            Rvec_n(:,ii)=sqrtm(dt*Dmtx)*randn(4,1);
        end
        
        
        n(2:5,:)=n(2:5,:)+dt*trans_n+Rvec_n;
        n(1,:)=ones(1,nsim)-sum(n(2:5,:),1);
        
        v=v-dt*Imemb;
    end
    fprintf('time %g ms\n',t)
end
realt=toc();
fprintf('ISIsHH85DAG realtime: %g sec\n',realt)

clf
plot(dt:dt:Tstop,vrec)

