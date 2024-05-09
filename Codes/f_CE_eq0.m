%Functions needed in turbulent boundary layer
function CE_eq0 = f_CE_eq0(nu, u, theta, H) 
    if(H<4) 
        H1 = 3.15 + 1.72/(H-1) - 0.01*(H-1)^2
    else
        H1 = 1.75 - 5.22273*H/(H+5.818181);
    end
    
    Re_t = abs(u*theta/nu)
    
    Cf0 = 0.01013/(log10(Re_t) - 1.02) - 0.00075;
    
    if(Cf0<0) 
        Cf0 = 1e-5;
    end
    H0 = 1/(1-6.55* sqrt(Cf0/2));
    
    Cf = Cf0*(0.9/(H/H0 - 0.4) - 0.5);
    
    tududx_eq0 = (1.25/H)*(Cf/2 - ((H-1)/(6.432*H))^2);
    CE_eq0 = H1*(Cf/2 - (H+1)*tududx_eq0);
    
end










