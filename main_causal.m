%% parameter initialization
K = 3;
A_all = binary(0:(2^K-1), K);
A_true = A_all;

% conditional distributions
n_in = size(A_true, 1);

G = [[1 0 0]; [1 1 0]; [1 0 1]; [1 1 1]; [1 0 1]; [0 1 0]; [0 1 1]; [0 0 0]];
[J, ~] = size(G);

theta_true = [
    [0.1 0.1 0.1 0.1 0.9 0.9 0.9 0.9];
    [0.05 0.05 0.4 0.4 0.7 0.7 0.95 0.95];
    [0.05 0.7 0.05 0.7 0.4 0.95 0.4 0.95];
    1-[0.1 1 2 3 4 5 6 6.9]/7;
    [-2 -0.5 -2 -0.5 2 0.5 2 0.5];          % Normal mean
    [-1.5 -1.5 1.5 1.5 -1.5 -1.5 1.5 1.5];  % Normal mean
    [0.5 -0.5 2 -2 0.5 -0.5 2 -2];          % Cauchy mean
    [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];      % Cauchy mean
];

theta_binary = [theta_true(1:4, :);
    1 - normcdf(0, theta_true(5, :), 1);
    1 - normcdf(0, theta_true(6, :), 1);
    0.5 + (1/pi) * atan(theta_true(7, :));
    0.5 + (1/pi) * atan(theta_true(8, :));
];

% latent proportion parameters
p1 = 0.66; p0 = 1-p1;
q00 = 0.8; q01 = 0.6; q10 = 0.4; q11 = 0.2;

% latent causal DAG Lambda (choose accordingly)
type = "collider";

if(type == "linear")
    nu_true = [p0*p1^2, p0^2*p1, p0^3, p0^2*p1, p1^2*p0, p1*p0^2, p1^2*p0, p1^3];
elseif(type == "collider")
    nu_true = [(1-q00)*p0^2, (1-q01)*p0*p1, q00*p0^2, q01*p0*p1, (1-q10)*p0*p1, ...
        (1-q11)*p1^2, q10*p0*p1, q11*p1^2];
elseif(type == "dependent")
    nu_true = [p0*p1^3+p1*p0^3, ones(1,3)*(2*p0^2*p1^2), ones(1,3)*(p1*p0^3+p0*p1^3), p0^4+p1^4];
end



%% penalized EM
K = 3;
M = 2^K;

lambda3_list = [1,10,100,1000];
tau2_list = [0.05 0.1];
rho2_list = 1;

N_vec = [1000 5000 10000];
epsilon = 0.125;
Nrep = 100;
Acc_G = zeros(Nrep, length(N_vec));
time_est = zeros(Nrep, length(N_vec));

for n_ind = 1:length(N_vec)
    N = N_vec(n_ind);

    pi_list = cell(Nrep,1);
    theta_list = cell(Nrep,1);
    parameter_list = cell(Nrep,1);
    
    parfor(n = 1:Nrep, 4)
        rng(n);
        [X_0, ~] = generate_X(N, nu_true, G, theta_true, A_true);
        X = X_0 > mean(X_0); % > 0;
        best_bic = Inf;
    
        tic;
        for j = 1:length(lambda3_list)
            for k = 1:length(tau2_list)
                for l = 1:length(rho2_list)
                    lambda3 = lambda3_list(j);
                    rho2 = rho2_list(l);
                    tau2 = tau2_list(k);
                    ite = 1;
                    e = 10;
    
                    theta = rand(J,M)*0.3 + 0.7*theta_binary;
                    d = get_d(theta);
                    pi = 0.3*drchrnd(3*ones(1,M),1) + 0.7*nu_true;
                    eps = 5*1e-3;
                    MAXITERS = 10;
    
                    while e > eps && ite < MAXITERS
                        pi1 = pi;
                        [s,pi] = update_s_pi(X, pi, theta, 0);
                        [theta, d] = ADMM(X, d, theta, s, lambda3, rho2, tau2);
                        e = norm(pi1 - pi);
                        ite = ite + 1;
                    end          
    
                    bic = BIC(X,pi,theta);
    
                    if bic <= best_bic
                        pi2 = pi;
                        theta2 = theta;
                        s2 = s;
                        best_bic = bic;
                        best_lambda3 = lambda3;
                        best_rho2 = rho2;
                        best_tau2 = tau2;
                    end
                end
            end
        end
    
        pi = pi2;
        theta = theta2;
        s = s2;
    
        parameter = [best_lambda3, best_rho2, best_tau2];
        time_est(n) = toc;
    
        pi_list{n} = pi;
        theta_list{n} = theta;
        parameter_list{n} = parameter;
    end
    
    % estimate bipartite graph
    A = zeros(M,M);
    
    epsilon = 0.125;
    for n = 1:Nrep
        theta = theta_list{n};
        G_est = zeros(J, K);
        for j = 1:J
            for k = 1:K
                ind_1 = find(A_all(:,k) == 1);                             
                ind_0 = find(A_all(:,k) == 0);
                if median(abs(theta(j, ind_1) - theta(j, ind_0))) > epsilon
                    G_est(j,k) = 1;
                end
            end
        end
        Acc_G(n, n_ind) = sum(G_est ~= G, 'all');
    end

    % save latent proportions to estimate the latent DAG
    pii = zeros(Nrep, M);
    for i = 1:Nrep
        pi_tmp = pi_list{i};
        pii(i,:) = pi_tmp(1:max(M,end));
    end
    writematrix(pii, "pi_"+type+"_"+N+".csv")
    
    save("K=3_"+type+"_"+N+".mat");
end

% accuracy of the bipartite graph Gamma
mean(Acc_G)
