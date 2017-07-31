%The example of using:
%[x y] = shooting_method(@fce,0.001,1e-6,0,1,[1 3],'fd')  
%It is meant that will be solved the BVP ODE described in the function fun,
%on the interval (0,1) with boundary conditions y(0) = 1 and y'(1) = 3

function [eta, v,eigval] = shooting_method(fun,h,zero,a,b,con,type,init)  
    tic;
    if (h <= 0) || (zero <= 0) || (a>=b) || (length(con) ~= 2) || (length(type) ~= 2)
        error('Check validity of input parameters - non-negative of accuracy and length of arrays');
    end

    if nargin == 8
        shoot1 = init(1); shoot2 = init(2);
    else
        shoot1 = -10; shoot2 = 10;
    end
    
    if (type(1)=='f')
        a1 = [con(1) shoot1];
        a2 = [con(1) shoot2];
    else
        a1 = [shoot1 con(1)];
        a2 = [shoot2 con(1)];
    end  

    cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/shooting_gotler'
    % sovle for base flow 
    C=0.509; Pr=1; D=1; 
    % const is Gstar-Q
    const=1; deta = h;
    [eta,baseT,baseTdash,~,~] = baseflow(C,Pr,D,h,a,b);
    baseT = interp1(eta,baseT,a:h:b,'spline');
    baseTdash = interp1(eta,baseTdash,a:h:b,'spline');
    eta=a:h:b;
    % input an initial condition
    khat=1;
    kappa=khat;
    const=1;
    [eta, v,eigval] = shooting_gotler3(@gotler,h,a,b,khat);
    % calculate beta from shooting method
    beta=eigval;

cd '/Users/samtomlinson/Documents/CDT_year_1/MRESproject/Codes/evolution_system/Rayleigh_IC'

%     size(baseT)
%     size(baseTdash)
%     [~, F1] = RungeKutta(a,b,h,a1,fun,kappa,beta,const,baseT,baseTdash); 
%     [~, F2] = RungeKutta(a,b,h,a2,fun,kappa,beta,const,baseT,baseTdash);         
%     
%     if (type(2)=='f')
%         F1 = F1(1,end) - con(2);
%         F2 = F2(1,end) - con(2);
%         r = 1;
%     else
%         F1 = F1(2,end) - con(2); 
%         F2 = F2(2,end) - con(2);
%         r = 2;
%     end    
%     
%     if (F1*F2 > 0) 
%         error('The root of F function does not exist, for selected initialization parameters. Please, change the init array.')
%     end
%     
%     F3 = F1;
%     
%     while (abs(F3) > zero) 
%         
%         shoot3 = (shoot1 + shoot2)/2; 
%         
%         if (type(1)=='f')
%            a3 = [con(1) shoot3];            
%         else
%            a3 = [shoot3 con(1)];            
%         end           
%         
%         [x, F3] = RungeKutta(a,b,h,a3,fun,kappa,beta,const,baseT,baseTdash); 
%         y = F3; F3 = F3(r,end) - con(2); 
%         if (F1*F3 < 0)
%             shoot2 = shoot3; F2 = F3;            
%         elseif (F1*F2 < 0)
%             shoot1 = shoot3; F1 = F3;
%         else
%             error('Error');           
%         end
%     end  
%     
%     eigval = shoot3;
        
%     h = plot(x,y(1,:),'k-'); set(h,'linewidth',2);
%     hold on;
%     h = plot(x,y(2,:),'r-'); set(h,'linewidth',2);  
%     xlabel('{\it x}','FontSize',12);
%     ylabel('y({\it x }), y^{(1)}({\it x })','FontSize',12);
%     title('Solution of 1D Boundary Value Problem by Shooting Method','FontSize',12);    
%     set(gca,'FontSize',12);          
%     legend('Function','{1^{st}} Derivative','Location','Best');       
%     hold off;
    toc;