function [transition_upper, transition_lower] = transition(Vt, theta, X_airfoil, Y_airfoil, nu, num_upper_panel, num_lower_panel, stagnation) 
    %Divide upper lower for some parameters
    [XB_Upper, XB_Lower] = Divide_Upper_Lower(X_airfoil, stagnation); %we use the midpoint here
    [YB_Upper, YB_Lower] = Divide_Upper_Lower(Y_airfoil, stagnation); %we use the midpoint here
    [Vt_Upper, Vt_Lower] = Divide_Upper_Lower(Vt, stagnation);
    [theta_Upper, theta_Lower] = Divide_Upper_Lower(theta, stagnation);
    
    %Rotate for lower part of the data
    XB_Lower = RotateData(XB_Lower);
    YB_Lower = RotateData(YB_Lower);
    Vt_Lower = RotateData(Vt_Lower);
    theta_Lower = RotateData(theta_Lower);
    
   
    
    %calculating distance 
    xu = zeros(num_upper_panel,1); xl = zeros(num_lower_panel,1);
    
    %Upper
    for i = 1:num_upper_panel
        if(i == 1) 
            xu(i) = 0;
        else
           d = sqrt( (XB_Upper(i) - XB_Upper(i-1))^2 + (YB_Upper(i) - YB_Upper(i-1))^2 );
           xu(i) = xu(i-1) + d;
        end
    end
    
    %Lower
    for i = 1:num_lower_panel
        if(i == 1) 
            xl(i) = 0;
        else
           d = sqrt( (XB_Lower(i) - XB_Lower(i-1))^2 + (YB_Lower(i) - YB_Lower(i-1))^2 );
           xl(i) = xl(i-1) + d;
        end
    end
    
    
    %Michel method
    %upper : 
    transition_upper = num_upper_panel; upper = 0; i = 2;
    
    while(i<=num_upper_panel & upper == 0) 
        Re_theta = abs(Vt_Upper(i))*theta_Upper(i)/nu; %Reynolds number using momentum thickness
        Re_x = abs(Vt_Upper(i))*xu(i)/nu; %Reynolds number using position of point in body
        Re_transition = 2.9*(Re_x)^(0.4); %Transition Reynolds number using Michel method
        %Re_transition = 1.174*(1+22400/Re_x)*(Re_x)^0.46; %Transition Reynolds number using Cebeci-Smith transition
        
        if(Re_theta > Re_transition) 
            transition_upper = i;
            upper = 1;
        end
        
        i = i+1;
    end

    
    
    %lower :
    transition_lower = num_lower_panel; lower = 0; i = 2;
    while(i<=num_lower_panel & lower == 0) 
        Re_theta = abs(Vt_Lower(i))*theta_Lower(i)/nu; %Reynolds number using momentum thickness
        Re_x = abs(Vt_Lower(i))*xl(i)/nu; %Reynolds number using position of point in body
        Re_transition = 2.9*(Re_x)^(0.4); %Transition Reynolds number using Michel method
        %Re_transition = 1.174*(1+22400/Re_x)*(Re_x)^0.46; %Transition Reynolds number using Cebeci-Smith transition
        
        if(Re_theta > Re_transition) 
            transition_lower = i;
            lower = 1;
        end
        
        i = i+1;
    end
    
end