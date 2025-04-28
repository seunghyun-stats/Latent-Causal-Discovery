function v = val(a, X, s, theta, j, k, d, u, gamma)
% this function outputs the value of the objective function for updating
% theta in Algorithm 1
%
% @param a: current theta value
% @param X: data matrix
% @param s: posterior of latent classes
% @param theta: theta in the previous step
% @param j: jth item
% @param k: kth latent attribute
% @param d: d in Algorithm 1
% @param u: u in Algorithm 1 
% @param gamma: tuning parameter in TLP
% 
% @return v: value of the objective function

    N = size(X,1);
    M = size(d,2);
    
    v1 = -sum(s(:,k) .* X(:,j))/N * log(a);
    v2 = -sum(s(:,k) .* (1-X(:,j)))/N * log(1-a);
   
    v3 = 0;
    for l = 1:(k-1)
        v3 = v3 + (d(j,l,k) - (theta(j,l) - a) + u(j,l,k))^2;
    end
    
    v4 = 0;
    for l = (k+1):M
        v4 = v4 + (d(j,k,l) - (a - theta(j,l)) + u(j,k,l))^2;
    end
    
    v = v1 + v2 + (v3+v4)*gamma/2;
end