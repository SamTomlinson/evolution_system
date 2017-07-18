jumpu0 = function();

% Base flow
[eta,baseT,baseTdash,baseU,baseUdash] = baseflow(C,Pr,D,deta,-deta/2,Leta+deta/2);

% Import amplitude from fortran code
tildeA = ;

% Calculate kernal integral
Ku = ;

% Evalulate combined integral 
int = ;

% Calculate jump
jumpu0 = ((-16*pi*Uc*(beta^4))/(alpha*(c^3)))*int;



