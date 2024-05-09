function [H, theta, lambda, delta_star, tau_wall, Cf, separation_upper, separation_lower, num_upper_panel, num_lower_panel, stagnation, W, XP_Upper, YP_Upper, XP_Lower, YP_Lower] = Thwaites_Method(num_panel, rho, miu, midpoint, panel_length, beta, Vt, alpha, X_body, Y_body)
    
    %Upper and Lower Panel division due to stagnation
    [stagnation, num_upper_panel, num_lower_panel] = detect_stagnation(num_panel, Vt);

    %Divide Upper and Lower airfoil : X_panel, Y_panel, panel_length, beta and tangential velocity
    [XP_Upper, XP_Lower] = Divide_Upper_Lower(midpoint(:,1), stagnation);
    [YP_Upper, YP_Lower] = Divide_Upper_Lower(midpoint(:,2), stagnation);
    [XB_Upper, XB_Lower] = Divide_Upper_Lower(X_body, stagnation);
    [YB_Upper, YB_Lower] = Divide_Upper_Lower(Y_body, stagnation);
    [Panel_Length_Upper, Panel_Length_Lower] = Divide_Upper_Lower(panel_length, stagnation);
    [Vt_Upper, Vt_Lower] = Divide_Upper_Lower(Vt, stagnation);
    %[beta_Upper, beta_Lower] = Divide_Upper_Lower(beta, separation);
    
    %Rotate for lower part :
    XB_Lower = RotateData(XB_Lower);
    YB_Lower = RotateData(YB_Lower);
    XP_Lower = RotateData(XP_Lower);
    YP_Lower = RotateData(YP_Lower);
    Panel_Length_Lower = RotateData(Panel_Length_Lower);
    Vt_Lower = RotateData(Vt_Lower);
    
    %Thwaites Method
    [theta_upper, lambda_upper, separation_upper, delta_star_upper, tau_wall_upper, Cf_upper, num_of_separation_upper, W_upper, H_upper] = Calc_Thwaites_for_UpperLower(miu, rho, abs(Vt_Upper), alpha, XP_Upper, YP_Upper, Panel_Length_Upper, num_upper_panel); 
    [theta_lower, lambda_lower, separation_lower, delta_star_lower, tau_wall_lower, Cf_lower, num_of_separation_lower, W_lower, H_lower] = Calc_Thwaites_for_UpperLower(miu, rho, abs(Vt_Lower), alpha, XP_Lower, YP_Lower, Panel_Length_Lower, num_lower_panel);    
    
    lambda_lower = RotateData(lambda_lower);
    theta_lower = RotateData(theta_lower);
    Cf_lower = RotateData(Cf_lower);
    tau_wall_lower = RotateData(tau_wall_lower);
    delta_star_lower = RotateData(delta_star_lower);
    W_lower = RotateData(W_lower);
    H_lower = RotateData(H_lower);
    separation_lower = stagnation + 2 - separation_lower;
    
    
    %Combining data : 
    theta = zeros(num_panel,1); lambda = zeros(num_panel,1); delta_star = zeros(num_panel,1); tau_wall = zeros(num_panel,1); Cf = zeros(num_panel,1);
    W = zeros(num_panel,1); H = zeros(num_panel,1);
    
    for i = 1:num_panel
        if(i <=num_lower_panel)
            theta(i) = theta_lower(i);
            lambda(i) = lambda_lower(i);
            delta_star(i) = delta_star_lower(i);
            tau_wall(i) = tau_wall_lower(i);
            Cf(i) = Cf_lower(i);
            W(i) = W_lower(i);
            H(i) = H_lower(i);
        elseif (i>num_lower_panel)
            theta(i) = theta_upper(i - num_lower_panel);
            lambda(i) = lambda_upper(i - num_lower_panel);
            delta_star(i) = delta_star_upper(i - num_lower_panel);
            tau_wall(i) = tau_wall_upper(i - num_lower_panel);
            Cf(i) = Cf_upper(i - num_lower_panel);
            W(i) = W_upper(i - num_lower_panel);
            H(i) = H_upper(i - num_lower_panel);
        end
    end
    
    %First end
    theta(1) = theta(2);
    lambda(1) = lambda(2);
    delta_star(1) = delta_star(2);
    tau_wall(1) = tau_wall(2);
    Cf(1) = Cf(2);
    W(1) = W(2);
    H(1) = H(2);
    
    %last end
    theta(num_panel) = theta(num_panel-1);
    lambda(num_panel) = lambda(num_panel-1);
    delta_star(num_panel) = delta_star(num_panel-1);
    tau_wall(num_panel) = tau_wall(num_panel-1);
    Cf(num_panel) = Cf(num_panel-1);
    W(num_panel) = W(num_panel-1);
    H(num_panel) = H(num_panel-1);
    
    %at stagnation point
    stagnation
    theta(stagnation) = theta(stagnation+1);
    lambda(stagnation) = lambda(stagnation+1);
    delta_star(stagnation) = delta_star(stagnation-1);
    tau_wall(stagnation) = tau_wall(stagnation+1);
    Cf(stagnation) = Cf(stagnation+1);
    W(stagnation) = W(stagnation+1);
    H(stagnation) = H(stagnation+1);
    
    %separation location
    separation_upper = stagnation + separation_upper;
    
end