function [Quantity_Upper, Quantity_Lower] = Divide_Upper_Lower(Quantity, upper_lower_boundary)
    [n,u] = size(Quantity);
    
    lower = upper_lower_boundary ; upper = n - lower;
    Quantity_Upper = zeros(upper,1); Quantity_Lower = zeros(lower,1);
    
    for i = 1:n
        if(i<=upper_lower_boundary)
            Quantity_Lower(i) = Quantity(i);
        else
            Quantity_Upper(i-upper_lower_boundary) = Quantity(i);
        end
    end
end