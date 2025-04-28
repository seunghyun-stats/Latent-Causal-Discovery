function [theta, d] = ADMM(X, d, theta, s, lambda, gamma, tau)
% This function update the theta and d in ADMM
%
% @param X : data matrix
% @param d : difference between item parameter pairs
% @param theta: item parameter matrix
% @param s: posterior of latent classes
% @param lambda: penalty coefficient
% @param gamma: tuning parameter 
% @param tau: tuning parameter
%
% @return theta: updated item parameter matriix
% @return d: updated item parameter pair differences

    MAXITERS = 10;
    TOL = 1e-2;
    e1 = 1;
    ite = 1;
    
    [J,M] = size(theta);
    u = zeros(J,M,M);
    
    while e1 > TOL && ite < MAXITERS
        d0 = d;
        for j = 1:J
            for m = 1:M
                fun = @(x) val(x,X,s,theta,j,m,d,u,gamma);
                theta(j,m) = fminbnd(fun,0,1);
            end
        end
        
        [d, u] = update_d_u(d, u, theta, lambda, gamma, tau);
        e = (d0 - d).^2;
        e1 = sum(sum(sum(e)));
        ite = ite + 1;
    end

end