%3 equations for head-lag entrainment equation
function [dtheta_dx, dH_dx, dCE_dx] = f_3Eqs(nu, u, dudx, theta, H, CE)
    %disp(H)
    if(H<4) 
        H1 = 3.15 + 1.72/(H-1) - 0.01*(H-1)^2;
        dH_dH1 = -0.15*(H-1)^2/(1.72 + 0.02*(H-1)^3);
    else
        H1 = 1.75 - 5.22273*H/(H+5.818181);
        dH_dH1 =  0.15/(- 5.22273*5.818181/(H+5.818181)^2);
    end
    
    Re_t = u*theta/nu
    
    Cf0 = 0.01013/(log10(Re_t) - 1.02) - 0.00075;
    
    if(Cf0<0) 
        Cf0 = 1e-5;
    end
    H0 = 1/(1-6.55* sqrt(Cf0/2));
    
    Cf = Cf0*(0.9/(H/H0 - 0.4) - 0.5);
    
    tududx_eq0 = (1.25/H)*(Cf/2 - ((H-1)/(6.432*H))^2);
    CE_eq0 = H1*(Cf/2 - (H+1)*tududx_eq0);
    
    F = (0.02*CE + CE^2 + 0.8*Cf0/3)/(0.01+CE);
    Ct_eq0 = f_Ct(CE_eq0,Cf0);
    Ct  = f_Ct(CE,Cf0);
    
    Lambda = 1;
    
    C = Ct_eq0/(Lambda^2) - 0.32*Cf0;
    CE_eq = sqrt(C/1.2 + 0.0001) - 0.01;
    Ct_eq = f_Ct(CE_eq,Cf0);
    tududx_eq = (0.5*Cf - CE_eq/H1)/(H+1);
    
    %Output function
    dtheta_dx = 1*(0.5*Cf - (H+2)*theta*dudx/u);
    dH_dx = dH_dH1*(CE - H1*(0.5*Cf - (H+1)*theta*dudx/u))/theta;
    dCE_dx =  1*(F/theta)*((2.8/(H+H1))*(sqrt(Ct_eq)-Lambda*sqrt(Ct)) + tududx_eq - theta*dudx/u);
end





