% Finite differences code for solution of Gotler system 

%clear all; close all; clc;

%% create grid

Lx1=10; Leta=10;
Nx1 = 150; Neta = 100;
dx1 = Lx1/Nx1; deta = Leta/Neta;
x1 = (-Lx1/2-1.5*dx1:dx1:Lx1/2+0.5*dx1)'; eta = -deta/2:deta:Leta+deta/2;

%% change folder for gotler initial condtions

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/shooting_gotler'

% sovle for base flow 
C=0.509; Pr=1; D=1; 
% const is Gstar-Q
const=1;
[eta,baseT,baseTdash,baseU,baseUdash] = baseflow(C,Pr,D,deta,-deta/2,Leta+deta/2);

baseT = interp1(eta,baseT,-deta/2:deta:Leta+deta/2,'spline');
baseTdash = interp1(eta,baseTdash,-deta/2:deta:Leta+deta/2,'spline');
baseU = interp1(eta,baseU,-deta/2:deta:Leta+deta/2,'spline');
baseUdash = interp1(eta,baseUdash,-deta/2:deta:Leta+deta/2,'spline');
eta=-deta/2:deta:Leta+deta/2;

% input an initial condition
khat=1;
[eta, v1,eigval] = shooting_gotler3(@gotler,deta,-deta/2,Leta+deta/2,khat);
% calculate beta from shooting method
beta=eigval;

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/evolution_system'

%% change folder for rayeligh initial condtions

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/evolution_system/Rayleigh_IC'

% sovle for base flow 
C=0.509; Pr=1; D=1; 
% const is Gstar-Q
const=1;
[eta,baseT,baseTdash,baseU,baseUdash] = baseflow(C,Pr,D,deta,-deta/2,Leta+deta/2);

baseT = interp1(eta,baseT,-deta/2:deta:Leta+deta/2,'spline');
baseTdash = interp1(eta,baseTdash,-deta/2:deta:Leta+deta/2,'spline');
baseU = interp1(eta,baseU,-deta/2:deta:Leta+deta/2,'spline');
baseUdash = interp1(eta,baseUdash,-deta/2:deta:Leta+deta/2,'spline');
eta=-deta/2:deta:Leta+deta/2;

% input an initial condition
khat=1;
[eta, v2,eigval] = shooting_method(@fun,deta,0.01,-deta/2,Leta+deta/2,[0,0],'ff');
% calculate beta from shooting method
kappa=eigval;

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/evolution_system'

%% initial set-up

% preallocate solution vector
v0sol=zeros(length(eta),length(x1));
q0sol=zeros(length(eta),length(x1));
T0sol=zeros(length(eta),length(x1));
% u0sol=zeros(length(eta),length(x1));

% set up initial conditions to the left for gotler 
v0sol(:,1)=v1(1,:)*exp(beta*(-Lx1/2-1.5*dx1));
v0sol(:,2)=v1(1,:)*exp(beta*(-Lx1/2-0.5*dx1));
% u0sol(:,1)=-trapz(baseU.*v0sol(:,1)/baseT);
% u0sol(:,2)=-trapz(baseU.*v0sol(:,2)/baseT);
T0sol(:,1)=((-baseTdash.*v1(1,:))./(baseT.*khat))*exp(beta*(-Lx1/2-1.5*dx1));
T0sol(:,2)=((-baseTdash.*v1(1,:))./(baseT.*khat))*exp(beta*(-Lx1/2-1.5*dx1));
q0sol(:,1)=zeros(size(v1(1,:)));
q0sol(:,2)=(v0sol(:,2)-v0sol(:,1))/(dx1);

% set up initial conditions to the left for rayleigh aswell 

v0sol(:,1)=v0sol(:,1)+(v2(1,:)*exp(2*kappa*(-Lx1/2-1.5*dx1)))';
v0sol(:,2)=v0sol(:,2)+(v2(1,:)*exp(2*kappa*(-Lx1/2-0.5*dx1)))';
% u0sol(:,1)=-trapz(baseU.*v0sol(:,1)/baseT);
% u0sol(:,2)=-trapz(baseU.*v0sol(:,2)/baseT);
T0sol(:,1)=T0sol(:,1)+(((-baseTdash.*v2(1,:))./(baseT.*khat))*exp(beta*(-Lx1/2-1.5*dx1)))';
T0sol(:,2)=T0sol(:,2)+(((-baseTdash.*v2(1,:))./(baseT.*khat))*exp(beta*(-Lx1/2-1.5*dx1)))';
q0sol(:,1)=q0sol(:,1)+(zeros(size(v2(1,:))))';
q0sol(:,2)=q0sol(:,2)+(v0sol(:,2)-v0sol(:,1))/(dx1);


%% marching downstream

Nx=length(eta)-2;

for i = 2:Nx+1
    
    % march T0
    for j=1:Neta+2
        T0sol(j,i+1) = T0sol(j,i) - dx1*(baseTdash(j)/baseT(j))*v0sol(j,i);
    end
    
    % calculate matrix A
    b=1-2./((khat.^2).*(baseT.^2).*deta.^2);
    c=baseTdash./((khat.^2).*(baseT.^3).*2.*deta) ...
        + 1./((khat.^2).*(baseT.^2).*(deta^2));
    a=-baseTdash./((khat.^2).*(baseT.^3).*2.*deta) ...
        + 1./((khat.^2).*(baseT.^2).*(deta.^2));
    T = b.*diag(ones(length(eta),1)) + c.*diag(ones(length(eta)-1,1),1)...
        + a.*diag(ones(length(eta)-1,1),-1);
    
    % create matrix A
    A = zeros(Neta+2);
    for j = 2:Nx+1
    A(j, j-1) = -(1)/((beta^2)*(baseT(j)^2)*(deta^2)) ...
        - (2*baseTdash(j))/((beta^2)*(baseT(j)^3)*(2*deta));
    A(j, j) = 1 + (2)/((beta^2)*(baseT(j)^2)*(deta^2));
    A(j, j+1) =  -(1)/((beta^2)*(baseT(j)^2)*(deta^2)) ...
        + (2*baseTdash(j))/((beta^2)*(baseT(j)^3)*(2*deta));
    end
    % with bcs
    A(1,1)=0.5; A(1,2)=0.5; A(Nx+2,Nx+1)=0.5; A(Nx+2,Nx+2)=0.5;
    
    % create matrix B
    B = zeros(Neta+2);
    for j = 1:Nx+2
    B(j, j) = const*baseT(j)^(-1);
    end
    
    % Calculate q0
    q0sol(:,i+1)=(A\B)*T0sol(:,i+1);
    
    % Calculate v0
    v0sol(:,i+1) = v0sol(:,i) + dx1*q0sol(:,i);
    
    
end

contourf(x1,eta,v0sol)
colorbar

