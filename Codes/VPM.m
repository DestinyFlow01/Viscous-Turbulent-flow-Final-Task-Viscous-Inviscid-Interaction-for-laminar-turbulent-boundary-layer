function [Vt,Vn, Cp, gamma, lambda,b] = VPM(alpha, V_inf, W, I, J, K, L, A, num_panel, beta, phi, panel_length)
    %Create vector b and solve for the source and vortex
    b = zeros(num_panel+1,1);
    
    %Initialization of b :
    for i = 1:num_panel+1
        if(i ~= num_panel + 1)
            b(i) = 2*pi*(W(i) - V_inf*cos(beta(i)));
        else
           b(i) = -2*pi*V_inf*(sin(beta(1)) + sin(beta(num_panel))); 
        end
    end
    
    %Solving source and vortex strength
    x = inv(A)*b;
    
    %Source strength
    lambda = x(1:num_panel,1);
    %Vortex strength
    gamma = x(num_panel+1,1);
    
    
    %Calculate Vt Vn and Cp
    Vt = zeros(num_panel,1);
    Vn = zeros(num_panel,1);
    Cp = zeros(num_panel,1);
    
    %Calculating Velocity
    for i = 1:num_panel
        %Inisialisasi :
        Vn(i) = V_inf*cos(beta(i)) + lambda(i)/2;
        Vt(i) = V_inf*sin(beta(i)) + gamma/2;
        
        for j = 1:num_panel
            if(j~=i)
               Vn(i) = Vn(i) + (lambda(j)*I(i,j)-gamma*K(i,j))/(2*pi);
               Vt(i) = Vt(i) + (lambda(j)*J(i,j)-gamma*L(i,j))/(2*pi);
            end
        end
    end
    
    %Calculate pressure coefficient
    for i = 1:num_panel
        Cp(i) = 1-(Vn(i)^2+Vt(i)^2)/V_inf^2;
    end

end