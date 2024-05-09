function [stagnation, num_upper_panel, num_lower_panel] = detect_stagnation(num_panel, Vt)
    %detecting stagnation
    stagnation = 0;
    for i = 2:num_panel
        if(Vt(i)/Vt(i-1) <=0)
            stagnation = i - 1; %focus on lower part
            break
        end
    end
    
    %Dividing upper and lower region
    num_upper_panel = num_panel - stagnation;
    num_lower_panel = stagnation; 
    
end