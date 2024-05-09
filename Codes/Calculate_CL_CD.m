function [CL, CD, CD_Viscous, CD_Pressure] = Calculate_CL_CD(Cp, Cf, beta, phi, alpha, num_panel, panel_length)
    %normal and axial coefficient : 
    CN_Cp = zeros(num_panel,1);
    CN_Cf = zeros(num_panel,1);
    CA_Cp = zeros(num_panel,1);
    CA_Cf = zeros(num_panel,1);
    
    CL = 0; CD = 0; CD_Viscous = 0; CD_Pressure = 0; 
    
    for i = 1:num_panel
        %calculate from Cp
        CN_Cp(i) = -Cp(i)*panel_length(i)*sin(alpha + beta(i));
        CA_Cp(i) = -Cp(i)*panel_length(i)*cos(alpha + beta(i));
        
        %Calculate from Cf
        CA_Cf(i) = abs(Cf(i)*panel_length(i)*sin(beta(i)+alpha));
        CN_Cf(i) = CA_Cf(i)/tan(beta(i)+alpha);
        
        
        %Calculate CL, CD, CD_Viscous, CD_Pressure
        CL = CL + (CN_Cp(i) + CN_Cf(i))*cos(alpha) - (CA_Cp(i) + CA_Cf(i))*sin(alpha);
        CD_Viscous = CD_Viscous + CN_Cf(i)*sin(alpha) + CA_Cf(i)*cos(alpha);
        CD_Pressure = CD_Pressure + CN_Cp(i)*sin(alpha) + CA_Cp(i)*cos(alpha);
        CD = CD_Viscous + CD_Pressure;
        
    end
    
end