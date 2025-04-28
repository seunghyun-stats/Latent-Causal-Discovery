function [d, u] = update_d_u(d, u, theta, lambda, gamma, tau)
% this function updates d and u
% we only need to update them with i < j
%
% @param d: d in the previous iteration
% @param u: u in the previous iteration
% @param theta: item parameter matrix
% @param lambda: tuning parameter for TLP
% @param gamma: tuning parameter for TLP
% @param tau: threshold for TLP
%
% @return d: updated d
% @return u: updated u

    [J,M] = size(theta);
    
    for j = 1:J
        for k = 1:M
            for l = (k+1):M
                tmp = theta(j,k) - theta(j,l) - u(j,k,l);
                if abs(tmp) >= tau
                    d(j,k,l) = tmp;
                else
                    d(j,k,l) = ST(tmp, lambda/gamma);
                end
                u(j,k,l) = u(j,k,l) + d(j,k,l) - (theta(j,k) - theta(j,l));
            end
        end
    end
    
end