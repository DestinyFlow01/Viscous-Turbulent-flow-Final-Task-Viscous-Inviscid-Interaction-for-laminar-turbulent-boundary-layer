function [theta, lambda, num_panel_start_separate, delta_star, tau_wall, Cf, num_of_separation, W, H_lambda] = Calc_Thwaites_for_UpperLower(miu, rho, Vt, alpha, X_panel, Y_panel, panel_length, num_panel)
    %Calculation of Thwaites method main idea : 
    %From nu, U and body or panel point, calculate theta _sq and theta
    %Calculate lambda and look for separations at lambda < -0.09 and record the value
    %Calculate S(lambda) and H(lambda) and consequently tau_wall and delta_star
    %Further calculation is to calculate Cf for each panel 
    %Divide the condition for upper and lower
    
    %Transforming to body coordinate with AoA :
    X_panel_new = zeros(num_panel,1); Y_panel_new = zeros(num_panel,1);
    
    
    %num_panel is for upper or lower only
    for i = 1:num_panel
        X_panel_new(i) = X_panel(i)*cos(-alpha) - Y_panel(i)*sin(-alpha);
        Y_panel_new(i) = X_panel(i)*sin(-alpha) + Y_panel(i)*cos(-alpha);
    end
    
    %calculate theta_sq and theta
    nu = miu/rho;
    theta_sq = zeros(num_panel,1);
    for i = 1:num_panel-1
        integral = 0;
        for j = 1:i
            integral = integral + Vt(j)^5*abs(X_panel_new(j+1)-X_panel_new(j));
        end
        theta_sq(i) = 0.45*nu*integral/(Vt(i))^6;
    end
    
    theta = sqrt(theta_sq);
    
    %calculate lambda and find first separation point
    lambda = zeros(num_panel,1);
    num_of_separation = 0;
    num_panel_start_separate = 1;
    
    for i = 1:num_panel
        if(i ~= num_panel)
            lambda(i) = theta_sq(i)/nu*(Vt(i+1)-Vt(i))/(X_panel_new(i+1)-X_panel_new(i));
            
        elseif(i== num_panel)
            lambda(i) = lambda(i-1);
        end
        
        if(lambda(i) < -0.09)
            lambda(i) = -0.09;
            if(num_of_separation == 0)
                num_of_separation = num_of_separation+1;
                num_panel_start_separate = i
                X_panel_new(i)
            end
        end
    end
    
    %find separation point
    %{
    separation = zeros(num_of_separation,1);
    index_number = 0;
    
    for i = 1:num_panel
        if(lambda(i)<-0.09)
            index_number = index_number + 1;
            separation(index_number) = i;
        end
    end
    num_panel_start_separate = separation(1);
    %}
    
    %Calculate S(lambda) and H(lambda), tau_wall, delta_star and
    %transpiraton velocity and Cf
    S_lambda = zeros(num_panel,1);
    H_lambda = zeros(num_panel,1);  z = zeros(num_panel,1);
    tau_wall = zeros(num_panel,1);
    delta_star = zeros(num_panel,1);
    W = zeros(num_panel,1);
    Cf = zeros(num_panel,1);
   
    for i = 1:num_panel
        S_lambda(i) = (lambda(i) + 0.09)^(0.62);
        
        
        tau_wall(i) = miu*abs(Vt(i))*S_lambda(i)/theta(i);
        Cf(i) = tau_wall(i)/(0.5*rho*Vt(i)^2);
        
        
        
        z(i) = 0.25-lambda(i);
        H_lambda(i) = 2 + 4.14*z(i) - 83.5*z(i)^2 + 854*z(i)^3 - 3337*z(i)^4 + 4576*z(i)^5;
        
        if(H_lambda(i)<0)
            H_lambda(i) = 0;
        end
        
        delta_star(i) = theta(i)*H_lambda(i);
        W(i) = delta_star(i)*0.03;
        %{
        if(i ~= num_panel && X_panel_new(i+1)~=X_panel_new(i) )
            W(i) = (Vt(i+1)*delta_star(i+1) - Vt(i)*delta_star(i))/(X_panel_new(i+1)-X_panel_new(i));
        else
            W(i) = W(i-1);
        end
        %}
        
        
        
    end
    plot(H_lambda)
end