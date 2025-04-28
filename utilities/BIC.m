function v = BIC(X, pi, theta)
% This function outputs BIC 
%
% @param X : data matrix
% @param pi : proportion vector
% @param theta: item parameter matrix
%
% @return v : value of BIC

    N = size(X,1);
    K = length(pi);
    J = size(theta,1);
    
    % likelihood
    v1 = 0;
    for i = 1:N
        tmp = 0;
        for k = 1:K
            tmp = tmp + exp(log(pi(k)) + sum(X(i,:) .* log(theta(:,k)')) + sum((1 - X(i,:)) .* log(1 - theta(:,k)')));
            
        end
        v1 = v1 + log(tmp);
    end
    theta_round = round(theta,2);
    
    % penalty
    c = 0;
    for j = 1:J
        c = c + length(unique(theta_round(j,:)));
    end
    
    v2 = (K+c) * log(N);
    v = -v1 + v2;
end
