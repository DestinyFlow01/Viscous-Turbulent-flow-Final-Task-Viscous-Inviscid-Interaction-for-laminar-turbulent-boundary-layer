function Cf = f_Cf(nu, u, theta, H)
    Re_t = u*theta/nu
    
    Cf0 = 0.01013/(log10(Re_t) - 1.02) - 0.00075;
    
    if(Cf0<0) 
        Cf0 = 1e-5;
    end
    H0 = 1/(1-6.55* sqrt(Cf0/2));
    
    Cf = Cf0*(0.9/(H/H0 - 0.4) - 0.5);
    
end
