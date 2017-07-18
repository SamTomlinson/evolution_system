jumpT0 = function();

% Base flow
[eta,baseT,baseTdash,baseU,baseUdash] = baseflow(C,Pr,D,deta,-deta/2,Leta+deta/2);

% Import amplitude from fortran code
tildeA = ;

% Calculate kernal integral
Ku = ;

% Evalulate combined integral 
int = ;

% Calculate jump
jumpu0 = (2*cos(2*beta*z)*Re*pi*S1*Uc*Tc)/(alpha)*int;
