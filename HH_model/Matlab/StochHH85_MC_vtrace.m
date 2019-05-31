% Stochastic Hodgkin and Huxley model
% Voltage is shifted from original model to agree with current conventions (Vext = 0)
%
% This script will simulate 10 parallel simulations of 50 ms and plot the
% voltage traces
%
% Explicit Markov Chain calculations, using the Gillespie algorithm as
% employed by Chow&White 1996 (Biophys J)

nsim=10;
gK=36; gNa=120; gL=0.3; %mS/cm2
EK=-77; EL=-54.4; ENa=50; %mV
Tstop=50; dt=0.005; %ms
points = round(Tstop/dt);
NNa=3000;
NK=NNa*0.3;%NK=900;

threshold=-10;
Eff=[];

tic();

firetimes=[];
Na_trans=[-1 1 0 0 0 0 0 0
            1 -1 0 0 0 0 0 0
            0 -1 1 0 0 0 0 0
            0 1 -1 0 0 0 0 0
            0 0 -1 1 0 0 0 0
            0 0 1 -1 0 0 0 0
            -1 0 0 0 1 0 0 0
            1 0 0 0 -1 0 0 0
            0 -1 0 0 0 1 0 0
            0 1 0 0 0 -1 0 0
            0 0 -1 0 0 0 1 0
            0 0 1 0 0 0 -1 0
            0 0 0 -1 0 0 0 1
            0 0 0 1 0 0 0 -1
            0 0 0 0 -1 1 0 0
            0 0 0 0 1 -1 0 0
            0 0 0 0 0 -1 1 0
            0 0 0 0 0 1 -1 0
            0 0 0 0 0 0 -1 1
            0 0 0 0 0 0 1 -1]';

K_trans=[-1 1 0 0 0
        1 -1 0 0 0
        0 -1 1 0 0
        0 1 -1 0 0
        0 0 -1 1 0
        0 0 1 -1 0
        0 0 0 -1 1
        0 0 0 1 -1]';
xx=zeros(1,nsim);

vrec = zeros(points,nsim);
p=1;

v = -65*ones(1,nsim);
alpha_n=0.01*(v+55)./(1-exp(-(v+55)/10));
beta_n=0.125*exp(-(v+65)/80);
alpha_m=0.1*(v+40)./(1-exp(-(v+40)/10));
beta_m=4*exp(-(v+65)/18);
beta_h=(1)./(1+exp(-(v+35)/10));
alpha_h=0.07*exp(-(v+65)/20);
M=alpha_m./beta_m;
H=alpha_h./beta_h;
Nastatesum=(1+H).*(1+M).^3;
Nastates=round(NNa*[ones(1,nsim);3*M;3*M.^2;M.^3;H;3*M.*H;3*M.^(2).*H;M.^(3).*H]./(ones(8,1)*Nastatesum));
if sum(Nastates,1)~=ones(1,nsim)*NNa
    Nastates(1,:)=NNa-sum(Nastates(2:8,:),1);
end

Narates=[3*alpha_m.*Nastates(1,:)
beta_m.*Nastates(2,:)
2*alpha_m.*Nastates(2,:)
2*beta_m.*Nastates(3,:)
alpha_m.*Nastates(3,:)
3*beta_m.*Nastates(4,:)
alpha_h.*Nastates(1,:)
beta_h.*Nastates(5,:)
alpha_h.*Nastates(2,:)
beta_h.*Nastates(6,:)
alpha_h.*Nastates(3,:)
beta_h.*Nastates(7,:)
alpha_h.*Nastates(4,:)
beta_h.*Nastates(8,:)
3*alpha_m.*Nastates(5,:)
beta_m.*Nastates(6,:)
2*alpha_m.*Nastates(6,:)
2*beta_m.*Nastates(7,:)
alpha_m.*Nastates(7,:)
3*beta_m.*Nastates(8,:)];
next_evNa=-log(rand(1,nsim))./sum(Narates,1);

N=alpha_n./beta_n;
Kstatesum=(1+N).^4;
Kstates=round(NK*[ones(1,nsim);4*N;6*N.^2;4*N.^3;N.^4]./(ones(5,1)*Kstatesum));
Krates=[4*alpha_n.*Kstates(1,:)
beta_n.*Kstates(2,:)
3*alpha_n.*Kstates(2,:)
2*beta_n.*Kstates(3,:)
2*alpha_n.*Kstates(3,:)
3*beta_n.*Kstates(4,:)
alpha_n.*Kstates(4,:)
4*beta_n.*Kstates(5,:)];
next_evK=-log(rand(1,nsim))./sum(Krates,1);

spikes=[];
firing=zeros(1,nsim);
tint = 5;

for tt=dt:tint:Tstop
    for t = tt:dt:tt+tint-dt
        vrec(p,:) = v;
        p=p+1;

        Imemb=gK*Kstates(5,:).*(v-EK)/NK+gNa*Nastates(8,:).*(v-ENa)/NNa+gL*(v-EL);

        v=v-dt*Imemb;
        while any(t>=next_evNa)
            ii=find(t>=next_evNa);
            dist=cumsum(Narates./(ones(20,1)*sum(Narates,1)),1);
            ev=rand(1,nsim);
            for a=ii
                xx(a)=find(ev(a)<dist(:,a),1);
            end
            Nastates(:,ii)=Nastates(:,ii)+Na_trans(:,xx(ii));

            prev_ev=next_evNa(ii);
            alpha_m=0.1*(v+40)./(1-exp(-(v+40)/10));
            beta_m=4*exp(-(v+65)/18);
            beta_h=(1)./(1+exp(-(v+35)/10));
            alpha_h=0.07*exp(-(v+65)/20);
            Narates=[3*alpha_m.*Nastates(1,:)
            beta_m.*Nastates(2,:)
            2*alpha_m.*Nastates(2,:)
            2*beta_m.*Nastates(3,:)
            alpha_m.*Nastates(3,:)
            3*beta_m.*Nastates(4,:)
            alpha_h.*Nastates(1,:)
            beta_h.*Nastates(5,:)
            alpha_h.*Nastates(2,:)
            beta_h.*Nastates(6,:)
            alpha_h.*Nastates(3,:)
            beta_h.*Nastates(7,:)
            alpha_h.*Nastates(4,:)
            beta_h.*Nastates(8,:)
            3*alpha_m.*Nastates(5,:)
            beta_m.*Nastates(6,:)
            2*alpha_m.*Nastates(6,:)
            2*beta_m.*Nastates(7,:)
            alpha_m.*Nastates(7,:)
            3*beta_m.*Nastates(8,:)];
            next_evNa(ii)=prev_ev-log(rand(1,length(ii)))./sum(Narates(:,ii),1);
            %			pause
        end

        while any(t>=next_evK)
            ii=find(t>=next_evK);
            dist=cumsum(Krates./(ones(8,1)*sum(Krates,1)),1);
            ev=rand(1,nsim);
            for a=ii
                xx(a)=find(ev(a)<dist(:,a),1);
            end
            Kstates(:,ii)=Kstates(:,ii)+K_trans(:,xx(ii));
            
            prev_ev=next_evK(ii);
            alpha_n=0.01*(v+55)./(1-exp(-(v+55)/10));
            beta_n=0.125*exp(-(v+65)/80);
            Krates=[4*alpha_n.*Kstates(1,:)
            beta_n.*Kstates(2,:)
            3*alpha_n.*Kstates(2,:)
            2*beta_n.*Kstates(3,:)
            2*alpha_n.*Kstates(3,:)
            3*beta_n.*Kstates(4,:)
            alpha_n.*Kstates(4,:)
            4*beta_n.*Kstates(5,:)];
            next_evK(ii)=prev_ev-log(rand(1,length(ii)))./sum(Krates(:,ii),1);
            %			printf('%g\t%g\t%g\t%g\tt=%g\n',next_evK,v,alpha_n,beta_n,t)

        end
    end
    fprintf('time %g ms\n',t);
end

realt=toc();
fprintf('ISIsHH85CW realtime: %g sec\n',realt)

clf
plot(dt:dt:Tstop,vrec)
