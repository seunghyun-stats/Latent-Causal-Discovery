function [X, X_A] = generate_X(N, nu_true, Q, theta_true, A)

[J, K] = size(Q);
% n_in = size(nu_true, 1);
n_in = max(size(nu_true));

% generate multinomial counts
counts = mnrnd(N, nu_true);
X_A = zeros(N, K);
n = 1;
% X_A stores empirical CDF for categories: 1,2,3,...,2^K
for a = 1:length(nu_true)
    X_A(n:(n+counts(a)-1), 1:K) = repmat(A(a, 1:K), counts(a), 1);
    n = n+counts(a);
end
% permute A, size N by 1
X_A = X_A(randperm(N), 1:K);

ind_XA = bin2ind(X_A);




%%%%%%%%%% generate GDINA parameters %%%%%%%%%%%
p_correct = zeros(N, J);
for i = 1:N
    p_correct(i,:) = theta_true(:,ind_XA(i)+1)';
end

p_sub_1 = p_correct(:, 1:4,:); p_sub_2 = p_correct(:, 5:6); p_sub_3 = p_correct(:, 7:8); 
X = zeros(N,J);
X(:, 1:4) = double(rand(size(p_sub_1)) < p_sub_1);
X(:, 5:6) = normrnd(p_sub_2, 1, N, 2);
X(:, 7:8)= p_sub_3 + tan(pi * (rand(N,2) - 0.5));
end