% Finite differences code for solution of Gotler system, evaluates 
% intial conditions for Gotler modes and a pair of instability waves
% then marches the solution downstream to get the evolution of 
% velocoity and temperature.

clear all; close all; clc;

%% create grid

% INPUT : bounds
x1a=-10; x1b=20;
etaa=1; etab=10;
% lengths
Lx1=(x1b-x1a); Leta=etab-etaa;
% INPUT : number of grid points
Nx1 = 1000; Neta = 1000;
% stepsize
dx1 = Lx1/Nx1; deta = Leta/Neta;
% grid
x1 = (x1a-1.5*dx1:dx1:x1b+0.5*dx1)'; eta = etaa-deta/2:deta:etab+deta/2;


%% gotler initial condtions

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/shooting_gotler'

% paramters required in base flow
C=0.509; Pr=1; D=1; const=2;
% INPUT : spanwise wavenumber
khat=1;
% calculate solution using shooting method
[~, vg,eigval] = shooting_gotler3(@gotler,deta,etaa-deta/2,etab+deta/2,khat);
% calculate beta using shooting method
beta=sqrt(eigval);
% for checking whether eigenmodes look right
%plot(v1(1,:),eta)

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/evolution_system'

%% rayleigh initial condtions

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/shooting_gotler'

% paramters required in base flow
C=0.509; Pr=1; D=1; 
% INPUT : spanwise wavenumber
khat=1;
% calculate solution using shooting method
[~, vr,eigval] = shooting_gotler3(@gotler,deta,etaa-deta/2,etab+deta/2,khat);
% calculate beta using shooting method
kappa=sqrt(eigval)/4;
% for checking whether eigenmodes look right
%plot(v1(1,:),eta)

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/evolution_system'

%% base flow

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/shooting_gotler'

% calculate andresize base flow vectors
[eta,baseT,baseTdash,~,~]=baseflow(C,Pr,D,deta,etaa-deta/2,etab+deta/2);
baseT = interp1(eta,baseT,etaa-deta/2:deta:etab+deta/2,'spline');
baseTdash = interp1(eta,baseTdash,etaa-deta/2:deta:etab+deta/2,'spline');
eta = interp1(eta,eta,etaa-deta/2:deta:etab+deta/2,'spline');

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/evolution_system'

%% initial set-up

% preallocate solution vectors
v0sol=zeros(length(eta),length(x1));
q0sol=zeros(length(eta),length(x1));
T0sol=zeros(length(eta),length(x1));
% set up initial conditions to the for gotler 
v0sol(:,1)=vg(1,:)*exp(beta*(x1a-1.5*dx1));
v0sol(:,2)=vg(1,:)*exp(beta*(x1a-0.5*dx1));
T0sol(:,1)=((-baseTdash.*vg(1,:))./(baseT.*khat))*exp(beta*(x1a-1.5*dx1));
T0sol(:,2)=((-baseTdash.*vg(1,:))./(baseT.*khat))*exp(beta*(x1a-0.5*dx1));
q0sol(:,1)=beta*vg(1,:)*exp(beta*(x1a-1.5*dx1));
q0sol(:,2)=beta*vg(1,:)*exp(beta*(x1a-0.5*dx1));
% add initial conditions for rayleigh
v0sol(:,1)=v0sol(:,1)+(vr(1,:)*exp(2*kappa*(x1a-1.5*dx1)))';
v0sol(:,2)=v0sol(:,2)+(vr(1,:)*exp(2*kappa*(x1a-0.5*dx1)))';
T0sol(:,1)=T0sol(:,1)+(((-baseTdash.*vr(1,:))./(baseT.*khat))*exp(2*kappa*(x1a-1.5*dx1)))';
T0sol(:,2)=T0sol(:,2)+(((-baseTdash.*vr(1,:))./(baseT.*khat))*exp(2*kappa*(x1a-0.5*dx1)))';
q0sol(:,1)=2*kappa*vg(1,:)*exp(2*kappa*(x1a-1.5*dx1));
q0sol(:,2)=2*kappa*vg(1,:)*exp(2*kappa*(x1a-0.5*dx1));


%% marching downstream

% actual size for loops neglecting ghost nodes
Nx=length(eta)-2;
for i = 2:Nx+2
    % march T0sol
    for j=1:Neta+2
        T0sol(j,i+1) = T0sol(j,i) - dx1*(baseTdash(j)/baseT(j))*v0sol(j,i);
    end
    % create matrix A (see notes)
    A = zeros(Neta+2);
    for j = 2:Nx+1
    A(j, j-1) = -(1)/((beta^2)*(baseT(j)^2)*(deta^2)) ...
        - (2*baseTdash(j))/((beta^2)*(baseT(j)^3)*(2*deta));
    A(j, j) = 1 + (2)/((beta^2)*(baseT(j)^2)*(deta^2));
    A(j, j+1) =  -(1)/((beta^2)*(baseT(j)^2)*(deta^2)) ...
        + (2*baseTdash(j))/((beta^2)*(baseT(j)^3)*(2*deta));
    end
    % with bcs 
    A(1,1)=1-2/((beta^2)*(baseT(1)^2)*(deta^2))...
        -6*baseTdash(1)/((beta^2)*(baseT(1)^3)*(2*deta));
    A(1,2)=5/((beta^2)*(baseT(1)^2)*(deta^2))...
        +8*baseTdash(1)/((beta^2)*(baseT(1)^3)*(2*deta));
    A(1,3)=-4/((beta^2)*(baseT(1)^2)*(deta^2))...
        -2*baseTdash(1)/((beta^2)*(baseT(1)^3)*(2*deta));
    A(1,4)=1/((beta^2)*(baseT(1)^2)*(deta^2));
    A(Nx+2,Nx-1)=1/((beta^2)*(baseT(end)^2)*(deta^2));
    A(Nx+2,Nx)=-4/((beta^2)*(baseT(end)^2)*(deta^2))...
        +2*baseTdash(end)/((beta^2)*(baseT(end)^3)*(2*deta));
    A(Nx+2,Nx+1)=5/((beta^2)*(baseT(end)^2)*(deta^2))...
        -8*baseTdash(end)/((beta^2)*(baseT(end)^3)*(2*deta));
    A(Nx+2,Nx+2)=1-2/((beta^2)*(baseT(end)^2)*(deta^2))...
        +6*baseTdash(end)/((beta^2)*(baseT(end)^3)*(2*deta));
    % create matrix B (see notes)
    B = zeros(Neta+2);
    for j = 1:Nx+2
    B(j, j) = const*baseT(j)^(-1);
    end
    % Calculate q0sol
    q0sol(:,i+1)=(A\B)*T0sol(:,i+1);
    % Calculate v0sol
    v0sol(:,i+1) = v0sol(:,i) + dx1*q0sol(:,i+1); 
end

%% plotting

% contour plot of velocity evolution on whole domain 
contourf(x1,eta,v0sol,20)
colorbar


