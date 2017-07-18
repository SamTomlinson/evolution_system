% Finite differences code for solution of Gotler system 
%clear all; close all; clc;

% create grid
Lx1=10; Leta=10;
Nx1 = 100; Neta = 100;
dx1 = Lx1/Nx1; deta = Leta/Neta;
x1 = (-dx1:dx1:Lx1+dx1/2)'; eta = -deta/2:deta:Leta+deta/2;

% change folder
cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/shooting_gotler'
% sovle for base flow 
C=0.509; Pr=1; D=1;
[eta,baseT,baseTdash,baseU,baseUdash] = baseflow(C,Pr,D,deta,-deta/2,Leta+deta/2);

baseT = interp1(eta,baseT,-deta/2:deta:Leta+deta/2,'spline');
baseTdash = interp1(eta,baseTdash,-deta/2:deta:Leta+deta/2,'spline');
baseU = interp1(eta,baseU,-deta/2:deta:Leta+deta/2,'spline');
baseUdash = interp1(eta,baseUdash,-deta/2:deta:Leta+deta/2,'spline');
eta=-deta/2:deta:Leta+deta/2;

% input an initial condition
khat=1;
[eta, v,eigval] = shooting_gotler3(@gotler,deta,-deta/2,Leta+deta/2,khat);

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/evolution_system'

% start marching downstream 

% preallocate solution vector
v0sol=zeros(length(eta),length(x1));
q0sol=zeros(length(eta),length(x1)-1);
T0sol=zeros(length(eta),length(x1)-1);
u0sol=zeros(length(eta),length(x1));


% set up initial conditions to the left
v0sol(:,1)=v(1,:);
v0sol(:,2)=v(1,:);
u0sol(:,1)=-trapz(baseU.*v0sol(:,1)/baseT);
u0sol(:,2)=-trapz(baseU.*v0sol(:,2)/baseT);
T0sol(:,1)=(-baseTdash.*v(1,:))./(baseT.*khat);
T0sol(:,2)=(-baseTdash.*v(1,:))./(baseT.*khat);

for i = 2:Nx1+1
    
    % Matrix T
    b=1-2./((khat.^2).*(baseT.^2).*deta.^2);
    c=baseTdash./((khat.^2).*(baseT.^3).*2.*deta) ...
        + 1./((khat.^2).*(baseT.^2).*(deta^2));
    a=-baseTdash./((khat.^2).*(baseT.^3).*2.*deta) ...
        + 1./((khat.^2).*(baseT.^2).*(deta.^2));
    T = b.*diag(ones(length(eta),1)) + c.*diag(ones(length(eta)-1,1),1)...
        + a.*diag(ones(length(eta)-1,1),-1);
    
    % Calculate q0
    q0sol(:,i)=inv(T)*T0sol(:,i);
    
    % Calculate v0
    v0sol(:,i) = v0sol(:,i-1) + dx1*q0sol(:,i);
    
    % Calculate u0
    u0sol(:,i) = -trapz(baseU.*v0sol(:,i)/baseT);
    
    % Match T0
    T0sol(:,i+1) = T0sol(:,i) + dx1*(-baseTdash/baseT)*v0sol(:,i);
    
end

[X,Y] = meshgrid(x1,eta);
figure
quiver(X,Y,u0sol,v0sol)


