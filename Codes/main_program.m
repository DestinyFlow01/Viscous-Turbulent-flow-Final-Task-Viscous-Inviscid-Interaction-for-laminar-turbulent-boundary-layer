%Viscous Inviscid interaction using Thwaites Method
%{
   
    Made By : Bryan
              13619042

    Viscous Inviscid Interaction (VII) : 
    Laminar : Thwaites Method
    Turbulent : Head Lag - Entrainment Method
    Transition : Michel transition model

%}
clc; close all; clear all
%Flight condition
V_inf = 1;               %Freestream velocity
aoa = [0];%[-2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12]';              %AOA (degree)
miu = 0.0000181206;   % Dynamic Viscosity at SL
rho = 1.225;          % Density at SL
nu = miu/rho;         % Kinematic Viscosity
Re = V_inf/nu            % Reynolds number
alpha = aoa*pi/180;

%reading data from csv and build airfoil
M = readmatrix('naca2412.csv');
%flipping the coordinates : 
X_airfoil = M(:,1); Y_airfoil = -M(:,2);

%{
figure(1) 
plot(X_airfoil, Y_airfoil)
hold on 
scatter(X_airfoil, Y_airfoil)
%}

%making panel
[midpoint, panel_length, phi, num_panel] = Make_Panel(X_airfoil, Y_airfoil,M);

%Calculation for each aoa
[n_aoa, b] = size(aoa);

%force coefficients
CL = zeros(n_aoa,1);
CD = zeros(n_aoa,1);
CD_Pressure = zeros(n_aoa,1);
CD_Viscous = zeros(n_aoa,1);

%Pressure coefficient :
Cp = zeros(num_panel,n_aoa);

%Friction coefficient : 
Cf = zeros(num_panel,n_aoa);

%Boundary layer and VPM properties : 
H = zeros(num_panel,n_aoa);
theta = zeros(num_panel,n_aoa);
lambda = zeros(num_panel,n_aoa);
delta_star = zeros(num_panel,n_aoa);
tau_wall = zeros(num_panel,n_aoa);
separation_upper = zeros(1,n_aoa);
separation_lower = zeros(1,n_aoa);
Vt = zeros(num_panel,n_aoa);
Vn = zeros(num_panel,n_aoa);
W = zeros(num_panel,n_aoa);
transition_upper = zeros(n_aoa,1);
transition_lower = zeros(n_aoa,1);


for i = 1:n_aoa
    %Calculate beta (angles between aoa and panel normal direction)
    [phi,beta] = Beta(phi,alpha(i),num_panel);

    %Calculate I,J,K,L,A
    [I,J,K,L,A] = ComputeMatA(num_panel,X_airfoil, Y_airfoil, midpoint, phi, panel_length);

    %Iteration
    W = zeros(num_panel,1);
    error = 100;
    total_iteration = 0;
    iteration = 0;
    delta_star_temp = zeros(num_panel,1);
    

    while(iteration <= total_iteration && error > 10^-4)
        delta_star_old_temp = delta_star_temp;
        [Vt_temp,Vn_temp, Cp_temp, gamma, lambda_source,b] = VPM(alpha(i), V_inf, W, I, J, K, L, A, num_panel, beta, phi, panel_length);
        
        
        %For laminar case (Thwaites method)
        %[H_temp, theta_temp, lambda_temp, delta_star_temp, tau_wall_temp, Cf_temp, separation_upper_temp, separation_lower_temp, num_upper_panel, num_lower_panel, stagnation, W_temp, XP_Upper, YP_Upper, XP_Lower, YP_Lower] = Thwaites_Method(num_panel, rho, miu, midpoint, panel_length, beta, Vt_temp, alpha(i), X_airfoil, Y_airfoil);        
        
        %Detect transition (Michel transition model)
        %[transition_upper(i,1), transition_lower(i,1)] = transition(Vt_temp, theta_temp, midpoint(:,1), midpoint(:,2), nu, num_upper_panel, num_lower_panel, stagnation);
        
        %Turbulent case (Head Lag - Entrainment Method)
        %In head method, we solve 3 equations to obtain some boundary layer quantities
        %[Cf_temp, theta_temp, delta_star_temp, tau_wall_temp, W_temp] = Head_Lag_Method(nu, transition_upper(i), transition_lower(i), num_panel, midpoint,Vt_temp, theta_temp, H_temp, delta_star_temp, W_temp, Cf_temp, tau_wall_temp, rho, V_inf, num_lower_panel, stagnation, num_upper_panel);
        
        iteration = iteration+1;
        %error = sum(sum(abs(delta_star_old_temp-delta_star_temp)/delta_star_temp));
        
        %Transfering data to array
        Cp(:,i) = Cp_temp;                              % the same as 
        %Cf(:,i) = Cf_temp;                              % obtained via turbulent
        %theta(:,i) = theta_temp;                        % obtained via turbulent
        %lambda(:,i) = lambda_temp;                      % -
        %delta_star(:,i) = delta_star_temp;              % obtain via H
        %tau_wall(:,i) = tau_wall_temp;                  % obtained via turbulent
        %separation_upper(:,i) = separation_upper_temp;  %
        %separation_lower(:,i) = separation_lower_temp;  %
        Vt(:,i) = Vt_temp;                              %
        Vn(:,i) = Vn_temp;                              %
        %W(:,i) = W_temp;                                % obtained via turbulent
        %H(:,i) = H_temp;
        %transition_upper(i,1) = transition_upper(i,1) + stagnation;
        %transition_lower(i,1) = stagnation - transition_lower(i,1);

    end

    %Calculate CL, CD, CD_Viscous, CD_Pressure
    %[CL(i), CD(i), CD_Viscous(i), CD_Pressure(i)] = Calculate_CL_CD(Cp(:,i), Cf(:,i), beta, phi, alpha(i), num_panel, panel_length);
end



%{
%Plot aerodynamic force coefficients
figure(2)
plot(aoa, CL)
title('AOA vs CL')
xlabel('AOA'); ylabel('CL');

figure(3)
plot(aoa, CD)
title('AOA vs CD')
xlabel('AOA'); ylabel('CD');

%Plot boundary layer properties 
%Plot pressure coefficient 
figure(4)
plot(midpoint(:,1), -Cp(1:end,1)) 
title('Negative Pressure Coefficient (-C_p) for NACA 0012')
xlabel('x/c'); ylabel('-Cp')
legend('Simulation');

%Plot friction coefficient
figure(5)
plot(midpoint(:,1), Cf(1:end,1)) 
title('Friction Coefficient (C_f) for NACA 0012')
xlabel('x/c'); ylabel('C_f')
legend('Simulation')

%Plot Velocity Distribution 
figure(6)
V = sqrt(Vt.^2 + Vn.^2)/V_inf;
plot(midpoint(:,1), V(1:end,1)) 
title('Velocity distribution (v/V) for NACA 0012')
xlabel('x/c'); ylabel('v/V')
legend('Simulation')

%Plot Shear stress distribution at wall
figure(7)
plot(midpoint(:,1), tau_wall(1:end,1)) 
title('Wall shear stress (\tau_{wall}) for NACA 0012')
xlabel('x/c'); ylabel('\tau_{wall}')
legend('Simulation')

%Plot momentum thickness
figure(8)
plot(midpoint(:,1), theta(1:end,1)) 
title('Momentum thickness (\theta) for NACA 0012')
xlabel('x/c'); ylabel('\theta')

%Plot displacement thickness
figure(9)
plot(midpoint(:,1), delta_star(1:end,1)) 
title('Displacement thickness (\delta*) for NACA 0012')
xlabel('x/c'); ylabel('\delta')

%Calculating boundary layer thickness
delta = zeros(num_panel,n_aoa);
for n = 1:n_aoa
    for i = 1:num_panel
       if(i<=stagnation)
           delta(i,n) = Y_airfoil(i,1) - delta_star(i,n);
       else
           delta(i,n) = Y_airfoil(i,1) + delta_star(i,n);
       end
    end
end

%for all AOA
for n = 1:n_aoa
    figure(9+n) 
    plot(X_airfoil(1:end-1,1),Y_airfoil(1:end-1,1),'-',X_airfoil(1:end-1,1),delta(:,1),'-');
    hold on

    %Coordinate for separation
    sep_upper = separation_upper(1,n); sep_lower = separation_lower(1,n); 

    plot(X_airfoil(sep_upper,1),delta(sep_upper,1),'k*-',X_airfoil(sep_lower,1),delta(sep_lower,1),'k*-');
    %legend('Airfoil NACA0012','Boundary Layer','Separation Point at Upper BL','Separation Point at Lower BL','Location','southeast')
    title('Boundary Layer at AOA = ', aoa(n))
    axis([0 1 -0.09 0.09]);
    xlabel('X')
    ylabel('Y')            
    grid on; 
    hold off;
end
%}