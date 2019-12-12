%% 
% Model script for CW1 of ECM3739 Autumn 2018
% vars V,n
% free pars I0
% f=rhs(x,I0) where x(1)=V, x(2)=n
clear
format compact
%% Set your personal parameter value for EL
% insert your student number xxxyyyzzz (usually starts with a 6)
% rename as spiking_model_xxxyyyzzz.m and call this in your scripts for
% each questions
SNumber=660031764;
if SNumber==0
   error(['Enter your 9 digit student number in your local copy of spiking_model.m   ',...
            'See instruction sheet or contact James Rankin at j.a.rankin@exeter.ac.uk if unsure.']); 
end
rng(SNumber)
VK=-91-4*rand;
disp('Your personal value of parameter EL:');
format longg
disp(VK)
format short

%% Fixed parameters
C=20;
phi=1/8;

VCa=120;
VL=-60;
gK=8;
gL=2;
gCa=4;

V1=-2.8; % 1.2
V2=26; % 18
V3=12; % 12 
V4=17.4; % 17.4

%% nonlinear functions of V
boltz=@(V,Vh,Vt)(0.5)*(1+tanh((V-Vh)/Vt)); 
minf=@(V)boltz(V,V1,V2);
winf=@(V)boltz(V,V3,V4); 
tauw=@(V)1./cosh((V-V3)/(2*V4));

%% Define the right hand side rhs for the model x(1)=V, x(2)=n
I=@(V,w)gK*w*(V-VK)+gCa*minf(V).*(V-VCa)+gL*(V-VL);
rhs=@(x,I0)[1/C*(I0-I(x(1,:),x(2,:)));...
              phi*(winf(x(1,:))-x(2,:))/tauw(x(1,:))];


