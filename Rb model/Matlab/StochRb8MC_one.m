% Rubinstein model of Ranvier Node's action potentials as in
% Rubinstein JT (1995) Biophys J 68: 779-785.
% Mino H, Rubinstein JT, White JA (2002) Ann Biomed Eng 30: 578-587.
% Bruce IC (2007) Ann Biomed Eng 35: 315-318; 

% Voltage is shifted so that resting voltage = 0 mV (the reversal of the leak)

% This script simulates 1000 sweeps per stimulus and calculates
% firing efficiency, mean firing time and firing time variance
% 100 (nsim) simultaneous sweeps, 10 times
% Lines 69, 109, 110, 162 and 163 are intended to make a voltage trace
% It is advisable to reduce nsim before uncommenting them
%
% Coupled activation particles, 8-state Na channel
% Markov Chain modeling (Chow&White algorithm)

nsim=10;
gNa=0.02569*265; %mS/cm2
ENa=144; %mV
Cm=0.0000714*265; Rm=1953.49/265;
Tstop=1; dt=0.001; %ms
Idel=0; Idur=0.1; %µA, ms, ms
threshold=80;
points = round(Tstop/dt);
NNa=1000;  %Number of sodium channels


% Each row of this matrix represent one of the 20 possible transitions
transitions=[-1 1 0 0 0 0 0 0
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

Iamp=5.8;

Eff=[];
meanFT=[];
varFT=[];
xx=zeros(1,nsim);
NaNs=[];


tic()

p=1;
% Uncomment this only with a small nsim
vrec = zeros(points,nsim);
v = 0*ones(1,nsim);
alpha_m=1.872*(v-25.41)./(1-exp((25.41-v)/6.06));
beta_m=3.973*(21-v)./(1-exp((v-21)/9.41));
M=alpha_m./beta_m;
beta_h=(22.57)./(1+exp((56-v)/12.5));
alpha_h=-0.549*(27.74+v)./(1-exp((v+27.74)/9.06));
H=alpha_h./beta_h;

%calculation of steady state at t=0
statesum=(1+H).*(1+M).^3;
states=round(NNa*[ones(1,nsim);3*M;3*M.^2;M.^3;H;3*M.*H;3*M.^(2).*H;M.^(3).*H]./(ones(8,1)*statesum));
if sum(states,1)~=ones(1,nsim)*NNa
    states(1,:)=NNa-sum(states(2:8,:),1);
end
rates=[3*alpha_m.*states(1,:)
    beta_m.*states(2,:)
    2*alpha_m.*states(2,:)
    2*beta_m.*states(3,:)
    alpha_m.*states(3,:)
    3*beta_m.*states(4,:)
    alpha_h.*states(1,:)
    beta_h.*states(5,:)
    alpha_h.*states(2,:)
    beta_h.*states(6,:)
    alpha_h.*states(3,:)
    beta_h.*states(7,:)
    alpha_h.*states(4,:)
    beta_h.*states(8,:)
    3*alpha_m.*states(5,:)
    beta_m.*states(6,:)
    2*alpha_m.*states(6,:)
    2*beta_m.*states(7,:)
    alpha_m.*states(7,:)
    3*beta_m.*states(8,:)];
%time of the next transition (one per each simulation)
next_ev=-log(1-rand(1,nsim))./sum(rates,1);
firetime=zeros(1,nsim);
firing=zeros(1,nsim);

for t = dt:dt:Tstop
    vrec(p,:) = v;  %Uncomment this only with a small nsim
    p=p+1;
    if any(~firing&v>=threshold) %detection of action potentials
        ind=find(v>=threshold&~firing);
        firetime(ind)=t;
        firing(ind)=1;
    end
    
    Iapp=Iamp*(t>Idel&t<(Idel+Idur));
    Imemb=-Iapp+gNa.*(v-ENa).*states(8,:)/NNa+v./Rm;
    
    while any(t>=next_ev)  %the use of while (instead of if) allows more than one transition per time step
        %ii contains the indices of the simulations in which a transition is going to occur
        ii=find(t>=next_ev);
        dist=cumsum(rates./(ones(20,1)*sum(rates,1)),1);
        ev=rand(1,nsim);
        for a=ii
            xx(a)=find(ev(a)<dist(:,a),1);
        end
        states(:,ii)=states(:,ii)+transitions(:,xx(ii));
        
        prev_ev=next_ev(ii);
        alpha_m=1.872*(v-25.41)./(1-exp((25.41-v)/6.06));
        beta_m=3.973*(21-v)./(1-exp((v-21)/9.41));
        beta_h=(22.57)./(1+exp((56-v)/12.5));
        alpha_h=-0.549*(27.74+v)./(1-exp((v+27.74)/9.06));
        
        rates=[3*alpha_m.*states(1,:)
            beta_m.*states(2,:)
            2*alpha_m.*states(2,:)
            2*beta_m.*states(3,:)
            alpha_m.*states(3,:)
            3*beta_m.*states(4,:)
            alpha_h.*states(1,:)
            beta_h.*states(5,:)
            alpha_h.*states(2,:)
            beta_h.*states(6,:)
            alpha_h.*states(3,:)
            beta_h.*states(7,:)
            alpha_h.*states(4,:)
            beta_h.*states(8,:)
            3*alpha_m.*states(5,:)
            beta_m.*states(6,:)
            2*alpha_m.*states(6,:)
            2*beta_m.*states(7,:)
            alpha_m.*states(7,:)
            3*beta_m.*states(8,:)];
        % The next event is calculated from the last event (not from the current time)
        next_ev(ii)=prev_ev-log(1-rand(1,length(ii)))./sum(rates(:,ii),1);
    end
    v=v-dt*Imemb/Cm;
end

clf  %Uncomment this with a small nsim
plot(dt:dt:Tstop,vrec)


realt=toc();
fprintf('time = %g\n',realt)
