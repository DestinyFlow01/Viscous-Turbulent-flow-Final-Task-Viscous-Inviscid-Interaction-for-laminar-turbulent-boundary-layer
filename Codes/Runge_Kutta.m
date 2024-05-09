function [theta_next, H_next, CE_next] = Runge_Kutta(nu, dx, theta, H, CE, u, dudx, u_next, dudx_next) 
    u_half = 0.5*(u + u_next);
    dudx_half = 0.5*(dudx+dudx_next);
    
    %Grouping the variables : 
    y = zeros(3,1);
    y(1,1) = theta; y(2,1) = H; y(3,1) = CE;
    
    %Calculating Runge-Kutta parameters :
    y1 = y;
    k1 = dx*f_3Eqs(nu,u,dudx,y1(1,1), y1(2,1), y1(3,1));
    
    y2 = y + 0.5*k1;
    k2 = dx*f_3Eqs(nu,u_half,dudx_half,y2(1,1), y2(2,1), y2(3,1));
    
    y3 = y + 0.5*k2;
    k3 = dx*f_3Eqs(nu,u_half,dudx_half,y3(1,1), y3(2,1), y3(3,1));
    
    y4 = y + k3;
    k4 = dx*f_3Eqs(nu,u_next,dudx_next,y4(1,1), y4(2,1), y4(3,1));
    
    %variables of theta, H and CE in the next node : 
    y_next = y + (k1 + 2*k2 + 2*k3 +k4)/6;
    
    %Regrouping the variables
    theta_next = y_next(1,1);
    H_next = y_next(2,1);
    CE_next = y_next(3,1);
end