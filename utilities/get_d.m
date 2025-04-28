function d = get_d(theta)
% this function outputs differences between item parameter pairs
% 
% @param theta: item parameter matrix
%
% @return d: differences between item parameter pairs

    [J,M] = size(theta);
    d = zeros(J,M,M);
    
    for j = 1:J
        for k = 1:M
            for l = (k+1):M
                d(j,k,l) = theta(j,k) - theta(j,l);
            end
        end
    end
    
end