 function [midpoint, panel_length, phi, num_panel] = Make_Panel(X_airfoil, Y_airfoil,M)
    [num_panel,b] = size(M);
    num_panel = num_panel - 1;
    
    
    midpoint = zeros(num_panel,b); %midpoint coordinate (x,y)
    panel_length = zeros(num_panel,1); %panel length
    phi = zeros(num_panel,1); %angle of each panel 

    for i = 1:num_panel
        midpoint(i,1) = (X_airfoil(i) + X_airfoil(i+1))/2;
        midpoint(i,2) = (Y_airfoil(i) + Y_airfoil(i+1))/2;
        X = X_airfoil(i+1) - X_airfoil(i); Y = Y_airfoil(i+1) - Y_airfoil(i);
        panel_length(i) = sqrt(X^2+Y^2);
        phi(i) = atan2(Y,X);
    end
    %{
    scatter(midpoint(:,1), midpoint(:,2))
    title('Airfoil')
    legend('airfoil plot','airfoil point','midpoint')
    hold off
    %}

end