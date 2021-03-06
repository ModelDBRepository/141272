% Stochastic Hodgkin and Huxley model
% Voltage is shifted from original model to agree with current conventions (Vext = 0)
%
% This script will simulate 10 parallel simulations of 50 ms and plot the
% voltage traces
%
% Diffusion approximation algorithm with UNcoupled activation subunits, as
% implemented by Orio and Soudry 2011 (Submitted to PLoS One)

nsim=10;
gK=36; gNa=120; gL=0.3; %mS/cm2
EK=-77; EL=-54.4; ENa=50; %mV
NNa=300;
NK=NNa*0.3;%NK=900;
Tstop=50; dt=0.001; %ms
points = round(Tstop/dt);
Iamp=5.5:0.5:15;
Idel=1; Idur=1; %�A, ms, ms
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
m=ones(1,nsim)./(1+beta_m./alpha_m);
h=ones(1,nsim)./(1+beta_h./alpha_h);
tint=5;

for tt=dt:tint:Tstop
    for t = tt:dt:tt+tint-dt
        vrec(p,:) = v;
        p=p+1;

        Imemb=gK*n.^(4).*(v-EK)+gNa*m.^(3).*h.*(v-ENa)+gL*(v-EL);

        alpha_n=0.01*(v+55)./(1-exp(-(v+55)/10));
        beta_n=0.125*exp(-(v+65)/80);
        alpha_m=0.1*(v+40)./(1-exp(-(v+40)/10));
        beta_m=4*exp(-(v+65)/18);
        beta_h=(1)./(1+exp(-(v+35)/10));
        alpha_h=0.07*exp(-(v+65)/20);
        SDn = sqrt(abs(alpha_n.*(1-n)+beta_n.*n)/(dt*NK*4));
        SDm = sqrt(abs(alpha_m.*(1-m)+beta_m.*m)/(dt*NNa*3));
        SDh = sqrt(abs(alpha_h.*(1-h)+beta_h.*h)/(dt*NNa));

        n=n+dt*(alpha_n.*(1-n)-beta_n.*n+randn(1,nsim).*SDn);
        m=m+dt*(alpha_m.*(1-m)-beta_m.*m+randn(1,nsim).*SDm);
        h=h+dt*(alpha_h.*(1-h)-beta_h.*h+randn(1,nsim).*SDh);

        v=v-dt*Imemb;

    end
    fprintf('time %g ms\n',t)
end
realt=toc();
fprintf('ISIsHH2F realtime: %g sec\n',realt)

clf
plot(dt:dt:Tstop,vrec)

