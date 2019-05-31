// Rubinstein model of Ranvier Node's action potentials as in
// Rubinstein JT (1995) Biophys J 68: 779-785.
// Mino H, Rubinstein JT, White JA (2002) Ann Biomed Eng 30: 578-587.
// Bruce IC (2007) Ann Biomed Eng 35: 315-318; 

// Voltage is shifted so that resting voltage = 0 mV (the reversal of the leak)

// This scripts runs the deterministic version of the model, both in its
// coupled (8-state channel) and uncoupled (2-state) particles approach
// to show that in  the deterministic limit they are equivalent

stacksize('max');
gNa=0.02569; //mS/cm2
//gNa=0
ENa=144; //mV
Cm=0.0000714; Rm=1953.49;
Tstop=1; dt=0.0001; //ms
//Iamp=10;
Iamp=6/265;
Idel=0; Idur=0.1; //µA, ms, ms
threshold=80
funcprot(0);
function y=alpha_m(v),y=1.872*(v-25.41)./(1-exp((25.41-v)/6.06)),endfunction;
function y=alpha_h(v),y=-0.549*(27.74+v)./(1-exp((v+27.74)/9.06)),endfunction;
function y=beta_m(v),y=3.973*(21-v)./(1-exp((v-21)/9.41)),endfunction;
function y=beta_h(v),y=(22.57)./(1+exp((56-v)/12.5)),endfunction;
funcprot(1);


nsim=2;
p=1;
points = round(Tstop/(dt*recp))
vrec = zeros(points,nsim);
vars=zeros(points,4);
v=zeros(1,nsim);
m=1/(1+beta_m(v(1))/alpha_m(v(1)));
h=1/(1+beta_h(v(1))/alpha_h(v(1)));
M2=alpha_m(v(2))/beta_m(v(2));
H2=alpha_h(v(2))/beta_h(v(2));
statesum=(1+H2)*(1+M2)^3;
states=[1;3*M2;3*M2^2;M2^3;H2;3*H2*M2;3*H2*M2^2;H2*M2^3]/statesum;
if sum(states)<>1 then
    printf("alerta!");
end

firetime=zeros(1,nsim);
firing=zeros(1,nsim);

for t = dt:dt:Tstop
    vrec(p,1)=v(1);
    vrec(p,2)=v(2);
    vars(p,:)=[m h h*m^3 states(8)];
    p=p+1;

    //    if or(~firing&v>=threshold) then
    //        ind=find(v>=threshold&~firing);
    //        firetime(ind)=t;
    //        firing(ind)=1;
    //    end

    Iapp=Iamp*(t>Idel&t<(Idel+Idur));
    Imemb1=-Iapp+gNa*m^(3)*h*(v(1)-ENa)+v(1)/Rm;
    Imemb2=-Iapp+gNa*states(8)*(v(2)-ENa)+v(2)/Rm;

    m=m+dt*(alpha_m(v(1)).*(1-m)-beta_m(v(1)).*m);
    h=h+dt*(alpha_h(v(1)).*(1-h)-beta_h(v(1)).*h);

    alpm=alpha_m(v(2));betm=beta_m(v(2));alph=alpha_h(v(2));beth=beta_h(v(2));
    transition=[-3*alpm-alph,betm,0,0,beth,0,0,0;
    3*alpm,-2*alpm-betm-alph,2*betm,0,0,beth,0,0;
    0,2*alpm,-2*betm-alpm-alph,3*betm,0,0,beth,0;
    0,0,alpm,-3*betm-alph,0,0,0,beth;
    alph,0,0,0,-beth-3*alpm,betm,0,0;
    0,alph,0,0,3*alpm,-betm-2*alpm-beth,2*betm,0;
    0,0,alph,0,0,2*alpm,-2*betm-alpm-beth,3*betm;
    0,0,0,alph,0,0,alpm,-3*betm-beth];
    states=states+dt*transition*states;
    states(1)=1-sum(states(2:8));

    v(1)=v(1)-dt*Imemb1/Cm;
    v(2)=v(2)-dt*Imemb2/Cm;

end

scf(0)
clf
plot(dt:dt:Tstop,vrec)
scf(1)
clf
plot(dt:dt:Tstop,vars)
legend("$m$","$h$","$m^{3}h$","$m_3h_1$")
firetime

