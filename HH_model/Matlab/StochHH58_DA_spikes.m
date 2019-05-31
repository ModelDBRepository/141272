% Stochastic Hodgkin and Huxley model
% Voltage is shifted from original model to agree with current conventions (Vext = 0)
%
% This script will simulate 500 seconds (10 parallel simulations of 50 s)
% and output a .cvs file with the spike intervals
%
% Diffusion approximation algorithm, as implemented by Orio and Soudry 2011
% (Submitted to PLoS One)

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

        ma=abs(m);
        R=randn(10,nsim).*sqrt((dt/NNa)*[3*am.*ma(1,:)+bm.*ma(2,:);
        2*am.*ma(2,:)+2*bm.*ma(3,:);
        am.*ma(3,:)+3*bm.*ma(4,:);
        ah.*ma(1,:)+bh.*ma(5,:);
        ah.*ma(2,:)+bh.*ma(6,:);
        ah.*ma(3,:)+bh.*ma(7,:);
        ah.*ma(4,:)+bh.*ma(8,:);
        3*am.*ma(5,:)+bm.*ma(6,:);
        2*am.*ma(6,:)+2*bm.*ma(7,:);
        am.*ma(7,:)+3*bm.*ma(8,:)]);
        Wtm=[-R(1,:)+R(2,:)+R(5,:);
        -R(2,:)+R(3,:)+R(6,:);
        -R(3,:)+R(7,:);
        R(8,:)-R(4,:);
        -R(8,:)+R(9,:)-R(5,:);
        -R(9,:)+R(10,:)-R(6,:);
        -R(10,:)-R(7,:)];
        m(2:8,:)=m(2:8,:)+dt*trans_m+Wtm;
        m(1,:)=ones(1,nsim)-sum(m(2:8,:),1);

        trans_n=[-(bn+3*an).*n(2,:)+4*an.*n(1,:)+2*bn.*n(3,:);
        -(2*bn+2*an).*n(3,:)+3*an.*n(2,:)+3*bn.*n(4,:);
        -(3*bn+an).*n(4,:)+2*an.*n(3,:)+4*bn.*n(5,:);
        -4*bn.*n(5,:)+an.*n(4,:)];

        na=abs(n);
        R=randn(4,nsim).*sqrt((dt/NK)*[4*an.*na(1,:)+bn.*na(2,:);
        3*an.*na(2,:)+2*bn.*na(3,:);
        2*an.*na(3,:)+3*bn.*na(4,:);
        an.*na(4,:)+4*bn.*na(5,:)]);
        Wtn=[-R(1,:)+R(2,:);
        -R(2,:)+R(3,:);
        -R(3,:)+R(4,:);
        -R(4,:);];
    
        n(2:5,:)=n(2:5,:)+dt*trans_n+Rvec_n;
        n(1,:)=ones(1,nsim)-sum(n(2:5,:),1);

        v=v-dt*Imemb;
    end
    fprintf('time %g ms\n',t)
end
realt=toc();
fprintf('ISIsHH85DA realtime: %g sec\n',realt)

ISI=[];
for a=1:nsim
    ISI=[ISI;diff(spikes(spikes(:,1)==a),2)];
end

csvwrite(sprintf('ISIsHH85DA-N%g-%gs.csv',NNa,realt),ISI);

