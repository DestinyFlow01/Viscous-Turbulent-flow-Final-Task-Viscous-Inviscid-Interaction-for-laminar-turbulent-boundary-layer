function [Data] = RotateData(Data)
    [a,b] = size(Data);
    Data_new = zeros(a,b);
    
    
    for i = 1:a
       Data_new(i) = Data(a - i + 1);
    end
    
    Data = Data_new;
end