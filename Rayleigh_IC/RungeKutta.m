function [x y] = RungeKutta(a,b,h,con,fun,kappa,beta,const,baseT,baseTdash)
    x = a:h:b; 
    n = length(x);
    y = zeros(length(con),n);
    y(:,1) = con;
    
    for k = 1:(n-1) 
        k1 = fun(x(k),y(:,k),kappa,beta,const,baseT(k),baseTdash(k));
        k2 = fun(x(k)+0.5*h,y(:,k)+0.5*h*k1',kappa,beta,const,baseT(k),baseTdash(k));
        k3 = fun(x(k)+0.5*h,y(:,k)+0.5*h*k2',kappa,beta,const,baseT(k),baseTdash(k));
        k4 = fun(x(k)+h,y(:,k)+h*k3',kappa,beta,const,baseT(k),baseTdash(k));
        y(:,k+1)=y(:,k)+h*(k1'+2*k2'+2*k3'+k4')/6;
    end