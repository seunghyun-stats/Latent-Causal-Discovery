function [s,pi] = update_s_pi(X,pi,theta,lambda)
% this function updates s and pi
%
% @param X: data matrix
% @param pi: proportion vector
% @param theta: item parameter matrix
% @param lambda: tuning parameter for log penalty
%
% @return s: updated posterior 
% @return pi: updated proportion vector

    N = size(X,1);
    % M = length(pi);
    
    pi = reshape(pi, 1, length(pi));
    s = pi .* exp(X*log(theta)+(1-X)*log(1-theta));
    s = s./sum(s,2);
    
    pi = sum(s,1)/N;
    % pi(pi<0) = 1e-5;
    % pi = pi/sum(pi);
    
end