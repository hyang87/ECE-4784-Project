%Corporate with Chuyao Feng
%constants
gk=36e-3;
gna=120e-3;
gl=0.3e-3;
Ek=-12;
Ena=115;
El=10.6;
Iin=0;
%Iin=5;  % For Problem 3
V0=-70;
Cm=1e-3;
AlphaM=0.1*((25-V0)/(exp((25-V0)/10)-1));
BetaM=4*exp(-V0/18);
AlphaN=0.01*((10-V0)/(exp((10-V0)/10)-1));
BetaN=0.125*exp(-V0/80);
AlphaH=0.07*exp(-V0/20);
BetaH=1/(exp((30-V0)/10)+1);

%Initial Condition
M0=AlphaM/(AlphaM+BetaM);
N0=AlphaN/(AlphaN+BetaN);
H0=AlphaH/(AlphaH+BetaH);
Ina=M0^3*gna*H0*(V0-Ena);
Ik=N0^4*gk*(V0-Ek);  
Il=gl*(V0-El);
Iion=Iin-Ik-Ina-Il;  
Vector=ones(100,2);
Vector2=ones(100,2);
V0=-70;

%Euler's method set up
t_in=0;
t_final=100;
steps=100*100;
dt=(t_final-t_in)/steps;
t(1)=t_in;
Vm(1)=V0;
M(1)=M0;
N(1)=N0;
H(1)=H0;

for n=1:steps
    
    %implement euler's method
    t(n+1)=t(n)+dt;
    Vm(n+1)=Vm(n)+dt*(Iion/Cm);
    M(n+1)=M(n)+dt*(AlphaM*(1-M(n))-BetaM*M(n));
    N(n+1)=N(n)+dt*(AlphaN*(1-N(n))-BetaN*N(n));
    H(n+1)=H(n)+dt*(AlphaH*(1-H(n))-BetaH*H(n));
    
    %update related variables
    AlphaM=0.1*((25-Vm(n+1))/(exp((25-Vm(n+1))/10)-1));
    BetaM=4*exp(-Vm(n+1)/18);
    AlphaN=0.01*((10-Vm(n+1))/(exp((10-Vm(n+1))/10)-1));
    BetaN=0.125*exp(-Vm(n+1)/80);
    AlphaH=0.07*exp(-Vm(n+1)/20);
    BetaH=1/(exp((30-Vm(n+1))/10)+1);
    Ina=M(n+1)^3*gna*H(n+1)*(Vm(n+1)-Ena);  
    Ik=N(n+1)^4*gk*(Vm(n+1)-Ek);       
    Il=gl*(Vm(n+1)-El); 
    
%     if mod(n,50) == 0   % for Problem 2
%         Iin = 5;       % for Problem 2
%     else               % for Problem 2
%         Iin =0;       % for Problem 2
%     end               % for Problem 2

    Iion=Iin-Ik-Ina-Il;      
end

%%Plots
figure;
plot(t,Vm);
figure;
plot(t,gna.*M.^3.*H,'r');
hold on
plot(t,gk.*N.^4,'b');



