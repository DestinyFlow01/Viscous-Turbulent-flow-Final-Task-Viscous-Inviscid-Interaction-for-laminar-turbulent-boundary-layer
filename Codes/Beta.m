function [phi,beta] = Beta(phi,alpha,num_panel)
    beta = zeros(num_panel,1);

    for i = 1:num_panel
       beta(i) = phi(i)-alpha + pi/2; 
       
       %making the angle to the range of 0 to 2*pi
       if(phi(i)<0)
          phi(i) = phi(i)+2*pi; 
       end
       
       if(beta(i)<0)
          beta(i) = beta(i)+2*pi; 
       end
       
    end
    
end