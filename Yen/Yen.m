% credit/license for this file: http://dismac.dii.unipg.it/paper/code.html from:
% Francesco Bianconi et al.: "A sequential machine vision procedure for
% assessing paper impurities" doi: 10.1016/j.compind.2013.12.001

function t = Yen(H)

%Computes image threshold based on Kapur's method%
%INPUT
%H        :   grey-scale histogram
%
%OUTPUT
%t        :   threshold

    n = length(H);
      
    %Grey values
    g = (0:n-1)';
    
    [w_b, w_f, mu_b, mu_f, var_b, var_f] = CumMeanVar(H);
    
    %Compute 'correlation'
    C_b = zeros(1,n);
    C_f = zeros(1,n);
    for i = 1:n
        
        H_b = [H(1:i); zeros(1,n-1)'];
        H_f = [zeros(1,i)'; H((i+1):n)];
        
        if w_b(i) ~= 0
            H_b = H/w_b(i);
        end
        if w_f(i) ~= 0
            H_f = H/w_f(i);
        end
        
        for j = 1:n
            if H_b(j) ~= 0
                C_b(i) = C_b(i) + H_b(j).* (-log2(H_b(j)));
            end
        end
        for j = 1:n
            if H_f(j) ~= 0
                C_f(i) = C_f(i) + H_f(j).* (-log2(H_f(j)));
            end
        end

    end
    
    [maxVal, indices] = max(C_b + C_f);
    t = g(indices(1));
    
end
