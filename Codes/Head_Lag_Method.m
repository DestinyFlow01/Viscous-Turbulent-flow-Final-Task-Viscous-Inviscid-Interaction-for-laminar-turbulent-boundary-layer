function [Cf, theta, delta_star, tau_wall, W] = Head_Lag_Method(nu, transition_upper, transition_lower, num_panel, midpoint,Vt, theta, H, delta_star, W, Cf, tau_wall, rho, V_inf, num_lower_panel, stagnation, num_upper_panel) 
    figure(1)
    plot(H)
    figure(2)
    plot(theta)

    X = midpoint(:,1); 
    for i = transition_upper + 1 : num_panel - 2
        dx = abs(X(i) - X(i-1));
        dudx1 = (Vt(i) - Vt(i-1))/abs(X(i) - X(i-1));
        dudx2 = (Vt(i+1) - Vt(i))/abs(X(i+1) - X(i));
        
        CE = f_CE_eq0(nu, Vt(i-1), theta(i-1), H(i-1))
        [theta(i), H(i), CE] = Runge_Kutta(nu, dx, theta(i), H(i), CE, Vt(i-1), dudx1, Vt(i), dudx2);
        
        Cf(i) = f_Cf(nu,Vt(i), theta(i), H(i));
        delta_star(i) = theta(i)*H(i);
    end
    
    for i = transition_lower : -1 : 1
        dx = abs(X(i) - X(i+1));
        dudx1 = (Vt(i) - Vt(i+1))/abs(X(i) - X(i+1));
        dudx2 = (Vt(i+1) - Vt(i))/abs(X(i+1) - X(i));
        
        CE = f_CE_eq0(nu, Vt(i+1), theta(i+1), H(i+1))
        [theta(i), H(i), CE] = Runge_Kutta(nu, dx, theta(i), H(i),CE, Vt(i+1), dudx1, Vt(i), dudx2);
        
        Cf(i) = f_Cf(nu,Vt(i), theta(i), H(i));
        delta_star(i) = theta(i)*H(i);
        
    end
    
    
end