function [I,J,K,L,A] = ComputeMatA(num_panel,X_airfoil, Y_airfoil, midpoint, phi, panel_length)
    %Translating variables to easier version
    X = X_airfoil; Y = Y_airfoil;
    x = midpoint(:,1); y = midpoint(:,2);
    
    
    %Creating matrix I,J,K,L
    I = zeros(num_panel, num_panel);
    J = I; K = I; L = I; 
    
    %Compute coefficients A,B,C,D,E
    for i = 1:num_panel-1
        for j = 1:num_panel-1
            if(j~=i) 
                A = -(x(i)-X(j))*cos(phi(j)) - (y(i)-Y(j))*sin(phi(j));
                B = (x(i)-X(j))^2 + (y(i)-Y(j))^2;
                Cn_S = sin(phi(i) - phi(j));
                Ct_S = -cos(phi(i) - phi(j));
                Dn_S = (y(i)-Y(j))*cos(phi(i)) - (x(i)-X(j))*sin(phi(i));
                Dt_S = (y(i)-Y(j))*sin(phi(i)) + (x(i)-X(j))*cos(phi(i));
                
                Cn_V = Ct_S;                                     
                Dn_V = Dt_S;   
                Ct_V = Cn_S;                                      
                Dt_V = -Dn_S;
                
                if (B<=A^2) 
                    E = 0;
                else
                    E = sqrt(B-A^2);
                end
                
                %compute I J K L
                if(E == 0 || E == NaN || E == Inf)
                    I(i,j) = 0;
                    J(i,j) = 0;
                    K(i,j) = 0;
                    L(i,j) = 0;
                else
                    %compute I (Source normal component)
                    term1 = 0.5*Cn_S*log((panel_length(j)^2 + 2*A*panel_length(j) + B)/B);
                    term2 = ((Dn_S-A*Cn_S)/E)*(atan2((A+panel_length(j)), E) - atan2(A,E));
                    I(i,j) = term1 + term2;
                    
                    %compute J (Source tangential component)
                    term1 = 0.5*Ct_S*log((panel_length(j)^2 + 2*A*panel_length(j) + B)/B);
                    term2 = ((Dt_S-A*Ct_S)/E)*(atan2((A+panel_length(j)),E) - atan2(A,E));
                    J(i,j) = term1 + term2;
                    
                    %compute K (Vortex normal component)
                    term1 = 0.5*Cn_V*log((panel_length(j)^2 + 2*A*panel_length(j) + B)/B);
                    term2 = ((Dn_V-A*Cn_V)/E)*(atan2((A+panel_length(j)),E) - atan2(A,E));
                    K(i,j) = term1 + term2;
                    
                    %compute L (Vortex tangential component)
                    term1 = 0.5*Ct_V*log((panel_length(j)^2 + 2*A*panel_length(j) + B)/B);
                    term2 = ((Dt_V-A*Ct_V)/E)*(atan2((A+panel_length(j)),E) - atan2(A,E));
                    L(i,j) = term1 + term2;
                    
                    %Problematic values :
                    if(imag(I(i,j)) ~= 0 || I(i,j) == Inf || I(i,j) == NaN)
                        I(i,j)=0;
                    end
                    if(imag(J(i,j)) ~= 0 || J(i,j) == Inf || J(i,j) == NaN)
                        J(i,j)=0;
                    end
                    if(imag(K(i,j)) ~= 0 || K(i,j) == Inf || K(i,j) == NaN)
                        K(i,j)=0;
                    end
                    if(imag(L(i,j)) ~= 0 || L(i,j) == Inf || L(i,j) == NaN)
                        L(i,j)=0;
                    end
                   
                end
            end
        end
    end
    
    %Create Matrix A
    A = zeros(num_panel + 1, num_panel + 1);
    
    for i = 1:num_panel+1
        for j = 1:num_panel+1
            if(i == j && i ~= num_panel + 1) 
                A(i,j) = pi;
            elseif(i == j && i == num_panel + 1)
                A(i,j) = 2*pi - (sum(L(1,:)) + sum(L(num_panel,:)));
            elseif(i~= num_panel+1 && j ~= num_panel+1 && i~=j)
                A(i,j) = I(i,j);
            elseif(i == num_panel+1)
                A(i,j) = J(1,j) + J(num_panel,j);
            elseif(j == num_panel + 1)
                A(i,j) = -sum(K(i,:));
            end
        end
    end
end