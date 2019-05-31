% Stochastic Hodgkin and Huxley model
% Voltage is shifted from original model to agree with current conventions (Vext = 0)
%
% This script will simulate 500 seconds (10 parallel simulations of 50 s)
% and output a .cvs file with the spike intervals
%
% Diffusion approximation algorithm, as described by Fox 1997 (Biophys J.)
% and implemented by Goldywn 2011 (Phys Rev E). The square root of  the
% diffusion matrix is calculated at each time step, using the steady state
% values of variables.

nsim=10;
gK=36; gNa=120; gL=0.3; %mS/cm2
EK=-77; EL=-54.4; ENa=50; %mV
Tstop=50000; dt=0.005; %ms
points = round(Tstop/dt);
NNa=3000;
NK=NNa*0.3;%NK=900;

threshold=-10;
Eff=[];

tic();

firetimes=[];
xx=zeros(1,nsim);

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
tint=50;
for tt=dt:tint:Tstop
    for t = tt:dt:tt+tint-dt
        
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
        
        %steady-state values of the variables
        mss=[bm.^(3).*bh;
            3*am.*bm.^(2).*bh;
            3*am.^(2).*bm.*bh;
            am.^(3).*bh;
            bm.^(3).*ah;
            3*am.*bm.^(2).*ah;
            3*am.^(2).*bm.*ah;
            am.^(3).*ah]./(ones(8,1)*((ah+bh).*(am+bm).^3));
        
        
        for ii = 1:nsim  %because this algorithm needs a matrix operation (line 90) it has to be done iterating over independent simulations
            ami=am(ii);bmi=bm(ii);ahi=ah(ii);bhi=bh(ii);
            Dmtx=[3*ami.*mss(1,ii)+(2*ami+bmi+ahi).*mss(2,ii)+2*bmi.*mss(3,ii)+bhi.*mss(6,ii),-2*(ami.*mss(2,ii)+bmi.*mss(3,ii)),0,0,-(ahi.*mss(2,ii)+bhi.*mss(6,ii)),0,0;
                -2*(ami.*mss(2,ii)+bmi.*mss(3,ii)),2*ami.*mss(2,ii)+(ami+2*bmi+ahi).*mss(3,ii)+3*bmi.*mss(4,ii)+bhi.*mss(7,ii),-(ami.*mss(3,ii)+3*bmi.*mss(4,ii)),0,0,-(ahi.*mss(3,ii)+bhi.*mss(7,ii)),0;
                0,-(ami.*mss(3,ii)+3*bmi.*mss(4,ii)),ami.*mss(3,ii)+(3*bmi+ahi).*mss(4,ii)+bhi.*mss(8,ii),0,0,0,-(ahi*mss(4,ii)+bhi.*mss(8,ii));
                0,0,0,ahi.*mss(1,ii)+(3*ami+bhi).*mss(5,ii)+bmi.*mss(6,ii),-(3*ami.*mss(5,ii)+bmi.*mss(6,ii)),0,0;
                -(ahi.*mss(2,ii)+bhi.*mss(6,ii)),0,0,-(3*ami.*mss(5,ii)+bmi.*mss(6,ii)),ahi.*mss(2,ii)+3*ami.*mss(5,ii)+(2*ami+bmi+bhi).*mss(6,ii)+2*bmi.*mss(7,ii),-2*(ami.*mss(6,ii)+bmi.*mss(7,ii)),0;
                0,-(ahi.*mss(3,ii)+bhi.*mss(7,ii)),0,0,-2*(ami.*mss(6,ii)+bmi.*mss(7,ii)),ahi.*mss(3,ii)+2*ami.*mss(6,ii)+(ami+2*bmi+bhi).*mss(7,ii)+3*bmi.*mss(8,ii),-(ami.*mss(7,ii)+3*bmi.*mss(8,ii));
                0,0,-(ahi*mss(4,ii)+bhi.*mss(8,ii)),0,0,-(ami.*mss(7,ii)+3*bmi.*mss(8,ii)),ahi.*mss(4,ii)+ami.*mss(7,ii)+(3*bmi+bhi).*mss(8,ii)]/NNa;
            
            Rvec(:,ii)=sqrtm(dt*Dmtx)*randn(7,1);
        end
        
        m(2:8,:)=m(2:8,:)+dt*trans_m+Rvec;
        m(1,:)=ones(1,nsim)-sum(m(2:8,:),1);
        
        
        trans_n=[-(bn+3*an).*n(2,:)+4*an.*n(1,:)+2*bn.*n(3,:);
            -(2*bn+2*an).*n(3,:)+3*an.*n(2,:)+3*bn.*n(4,:);
            -(3*bn+an).*n(4,:)+2*an.*n(3,:)+4*bn.*n(5,:);
            -4*bn.*n(5,:)+an.*n(4,:)];
        
        nss=[bn.^4;
            4*an.*bn.^3;
            6*an.^2 .*bn.^2;
            4*an.^3 .*bn;
            an.^4]./(ones(5,1)*((an + bn).^4));
        
        for ii = 1:nsim  %because this algorithm needs a matrix operation (line 90) it has to be done iterating over independent simulations
            ani=an(ii);bni=bn(ii);
            Dmtx = [4*ani*nss(1,ii)+(3*ani+bni)*nss(2,ii)+2*bni*nss(3,ii),-(3*ani*nss(2,ii)+2*bni*nss(3,ii)),0,0;
                -(3*ani*nss(2,ii)+2*bni*nss(3,ii)),3*ani*nss(2,ii)+2*(ani+bni)*nss(3,ii)+3*bni*nss(4,ii),-(2*ani*nss(3,ii)+3*bni*nss(4,ii)),0;
                0,-(2*ani*nss(3,ii)+3*bni*nss(4,ii)),2*ani*nss(3,ii)+(ani+3*bni)*nss(4,ii)+4*bni*nss(5,ii),-(ani*nss(4,ii)+4*bni*nss(5,ii));
                0,0,-(ani*nss(4,ii)+4*bni*nss(5,ii)),ani*nss(4,ii)+4*bni*nss(5,ii)]/NK;
            
            Rvec_n(:,ii)=sqrtm(dt*Dmtx)*randn(4,1);
        end
        
        n(2:5,:)=n(2:5,:)+dt*trans_n+Rvec_n;
        n(1,:)=ones(1,nsim)-sum(n(2:5,:),1);
        
        v=v-dt*Imemb;
    end
    fprintf('time %g ms\n',t)
end
realt=toc();
fprintf('ISIsHH85DAGss realtime: %g sec\n',realt)


ISI=[];
for a=1:nsim
    ISI=[ISI;diff(spikes(spikes(:,1)==a),2)];
end

csvwrite(sprintf('ISIsHH85DAGss-N%g-%gs.csv',NNa,realt),ISI);

