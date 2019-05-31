% Stochastic Hodgkin and Huxley model
% Voltage is shifted from original model to agree with current conventions (Vext = 0)
%
% This script will simulate 10 parallel simulations of 50 ms and plot the
% voltage traces
%
% Explicit Markov Chain with UNcoupled activation subunits, as
% described by Orio and Soudry 2011 (Submitted to PLoS One)

nsim=10;
gK=36; gNa=120; gL=0.3; %mS/cm2
EK=-77; EL=-54.4; ENa=50; %mV
Tstop=50; dt=0.001; %ms
points = round(Tstop/dt);
NNa=300;
NK=NNa*0.3;%NK=900;

Idel=1; Idur=1; %µA, ms, ms
threshold=-10;

vrec = zeros(points,nsim);
p=1;

tic();

v = -65*ones(1,nsim);
alpha_n=0.01*(v+55)./(1-exp(-(v+55)/10));
beta_n=0.125*exp(-(v+65)/80);
alpha_m=0.1*(v+40)./(1-exp(-(v+40)/10));
beta_m=4*exp(-(v+65)/18);
beta_h=(1)./(1+exp(-(v+35)/10));
alpha_h=0.07*exp(-(v+65)/20);        
n=ones(1,nsim)./(1+beta_n./alpha_n);
Nn=round(n*NK*4);
m=ones(1,nsim)./(1+beta_m./alpha_m);
Nm=round(m*NNa*3);
h=ones(1,nsim)./(1+beta_h./alpha_h);
Nh=round(h*NNa);
next_evn=-log(rand())./(Nn.*beta_n + (NK*4-Nn).*alpha_n);
next_evh=-log(rand())./(Nh.*beta_h + (NNa-Nh).*alpha_h);
next_evm=-log(rand())./(Nm.*beta_m + (NNa*3-Nm).*alpha_m);
tint=5;

for tt=dt:tint:Tstop
    for t = tt:dt:tt+tint-dt
        vrec(p,:) = v;
        p=p+1;
        
        m=Nm/(NNa*3); h=Nh/NNa; n=Nn/(NK*4);
        Imemb=gK*n.^(4).*(v-EK)+gNa*m.^(3).*h.*(v-ENa)+gL*(v-EL);

        v=v-dt*Imemb;

        while any(t>=next_evn)
            ii=find(t>=next_evn);
            Nn(ii)=Nn(ii)+1-2*(rand(1,length(ii))<=(Nn(ii).*beta_n(ii)./(Nn(ii).*beta_n(ii)+(NK*4-Nn(ii)).*alpha_n(ii))));
            prev_ev=next_evn(ii);
            alpha_n(ii)=0.01*(v(ii)+55)./(1-exp(-(v(ii)+55)/10));
            beta_n(ii)=0.125*exp(-(v(ii)+65)/80);
            next_evn(ii)=prev_ev-log(rand())./(Nn(ii).*beta_n(ii) + (NK*4-Nn(ii)).*alpha_n(ii));
        end

        while any(t>=next_evm)
            ii=find(t>=next_evm);
            Nm(ii)=Nm(ii)+1-2*(rand(1,length(ii))<=(Nm(ii).*beta_m(ii)./(Nm(ii).*beta_m(ii)+(NNa*3-Nm(ii)).*alpha_m(ii))));
            prev_ev=next_evm(ii);
            alpha_m(ii)=0.1*(v(ii)+40)./(1-exp(-(v(ii)+40)/10));
            beta_m(ii)=4*exp(-(v(ii)+65)/18);
            next_evm(ii)=prev_ev-log(rand())./(Nm(ii).*beta_m(ii) + (NNa*3-Nm(ii)).*alpha_m(ii));
        end

        while any(t>=next_evh)
            ii=find(t>=next_evh);
            Nh(ii)=Nh(ii)+1-2*(rand(1,length(ii))<=(Nh(ii).*beta_h(ii)./(Nh(ii).*beta_h(ii)+(NNa-Nh(ii)).*alpha_h(ii))));
            prev_ev=next_evh(ii);
            beta_h(ii)=(1)./(1+exp(-(v(ii)+35)/10));
            alpha_h(ii)=0.07*exp(-(v(ii)+65)/20);                       
            next_evh(ii)=prev_ev-log(rand())./(Nh(ii).*beta_h(ii) + (NNa-Nh(ii)).*alpha_h(ii));
        end    
    end
    fprintf('time %g ms\n',t)
end
realt=toc();
fprintf('ISIsHH2CW realtime: %g sec\n',realt)


clf
plot(dt:dt:Tstop,vrec)