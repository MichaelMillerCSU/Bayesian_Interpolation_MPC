clc
clear all
close all

% A = [1.2 2;
%      0 1.2];
% B = [0;
%      0.079];

% A = [0 1;
%      0.1 1.1];
% B = [0;
%      0.01];

%%%%%%%%%%%%%%%%%%
ra = randi(5000);
rng(3506) 
f_u =  @(t,x,u)(-[ -x(2,:) ; 1*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);
n = 2;
m = 1; % number of control inputs
deltaT = 0.05;
%Runge-Kutta 4
k1 = @(t,x,u) ( f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

%%%%%%%%%%%%%%%%%%


% A = [0.997561799256275	0.0995284960866619
% -0.0485855254061620	0.988369439858514];
% B = [0.00248409889321964
% 0.0489081095678455];

Nsim = 1000;
Ntraj = 100;

% liftFun = @(x) [x; x(1) + x(2); 2.*x(1) + x(2)];
liftFun = @(x) [x;];

Nlift =size( liftFun(zeros(2, 1)), 1);
Ubig = 20*rand([Nsim Ntraj]) - 10;
Xcurrent = (rand(2,Ntraj)*4 - 2);

[A, B, C, ~, ~, ~, ~] = VDP_System_Generate(liftFun, Ubig, Xcurrent, Nsim);

rank(ctrb(A, B))

[n, m] = size(B);

% A = [    A           B;
%      zeros(m, n) eye(m, m)]


%% Unconstrained MPC
Np = 20;
Nc = 20;


Z = [];
for i = 1 : Np
    Temp = [];
    for j = 1 : Np
        if j <= i
            Temp = [Temp A^(i - j) * B];
        else
            Temp = [Temp zeros(n, m)];
        end
    end
    Z = [Z; Temp];
end

A_bar = [];
for i = 1 : Np
    A_bar = [A_bar; A^i];
end

q = 1;
r = 1;

Q = q * eye(n, n);
R = r;
[K, P] = dlqr(A,B,Q,R);
Q_bar = kron(diag(ones(Np, 1)), Q);
Q_bar(end - n + 1 : end, end - n + 1 : end) = P;
R_bar = kron(diag(ones(Np, 1)), R);
x_init = [-0.55; 0.55];
% x_init = [1; -1];

x = liftFun(x_init);
f = 2 * x' * A_bar' * Q_bar * Z;
H = Z' * Q_bar * Z + R_bar;


state_independent_num = 100;
% control_independent_num = 30;

% N = control_independent_num * Nc * m;
% nn = state_independent_num * Np * n;

N = state_independent_num * Nc * m;
nn = state_independent_num * Np * n;

control_bound = 10;

Constrained_MPC = 1;
Unconstrained_MPC = 1 - Constrained_MPC;

if Np == Nc
    A_eq = [];
    b_eq = [];
else
    A_eq = zeros(Np - Nc, Np);
    b_eq = zeros(Np - Nc, 1);
    
    for i = 1 : Np - Nc
        A_eq(i, Nc + i - 1 : Nc + i) = [1 -1];
    end
end




%% Simulation 1
X_Set = [x];
U_Set = [];
J_Set = [];
J = 0;

for i = 1 : 200
    
    Q = q * eye(n, n);
    R = r;
    [K, P] = dlqr(A,B,Q,R);
    Q_bar = kron(diag(ones(Np, 1)), Q);
    Q_bar(end - n + 1 : end, end - n + 1 : end) = P;
    
    f = 2 * x' * A_bar' * Q_bar * Z;
    H = Z' * Q_bar * Z + R_bar;
    H = (H + H') / 2;
    if Constrained_MPC
%         disp('Constrained On')
        u_seq = quadprog(2 * H, f, [], [], A_eq, b_eq, -control_bound * ones(Np, 1), control_bound * ones(Np, 1));
    end
    
    if Unconstrained_MPC
%         disp('Unconstrained On')
        u_seq = -inv(2 * H) * f';
    end
    K_mpc = -inv(2 * H) * 2 * Z' * Q_bar * A_bar;
    K_mpc = K_mpc(1, :);
    abs(eig(A + B * K_mpc));
    w = 0.01 * rand(n, 1) - 0.005;
%     w = 0;
    
%     % Continue dynamics
%     load U_Save.mat;
%     u_seq = U_Set(:, i);
%     x_next = A * x + B * u_seq(1) + w;
    x_next = f_ud(0, C * x, u_seq(1));
%     J = J + x' * Q * x + u_seq(1)' * R * u_seq(1)';
    x = liftFun(x_next);
    X_Set = [X_Set x];
    U_Set = [U_Set u_seq(1)];
%     J_Set = [J_Set J];

end
    
figure(1)
plot(X_Set(1, :), 'LineStyle','--', 'LineWidth', 3.0) 
hold on


figure(2)
plot(X_Set(2, :), 'LineStyle','--', 'LineWidth', 3.0) 
hold on

% legend('$x_1$', '$x_2$', 'Interpreter', 'latex')


figure(3)
plot(U_Set(1, :), 'LineStyle','--', 'LineWidth', 3.0) 
hold on
% legend('$u$', 'Interpreter', 'latex')



% figure(4)
% plot(J_Set(1, :), 'LineStyle','-', 'LineWidth', 3.0) 
% hold on
% legend('$J$', 'Interpreter', 'latex')


figure (10)
plot(X_Set(1, :), 'LineStyle','--', 'LineWidth', 3.0) 
hold on
plot(X_Set(2, :), 'LineStyle','--', 'LineWidth', 3.0) 
% legend('$x_1$ without Bayesian', '$x_2$ without Bayesian', 'Interpreter', 'latex')


%% Simulation 2

%% Weights Offline Update
x = liftFun(x_init);
Weights_Update_On = 1;
r_set = [r];
q_set = [q];
q_r_relative_Set = [];
% XX = kron(ones(1, 50), x );
XX = 2 * rand(2, state_independent_num) - 1;
if Weights_Update_On
    iterative_step = 50;
    for j = 1 : iterative_step
        X_seq = [];
        U_seq = [];
        for k = 1 : state_independent_num
            Q = q * eye(n, n);
            R = r;
            [K, P] = dlqr(A,B,Q,R);
            Q_bar = kron(diag(ones(Np, 1)), Q);
%             Q_bar(end - n + 1 : end, end - n + 1 : end) = P;
%             Q_bar = kron(diag(ones(Np, 1)), Q);
            R_bar = kron(diag(ones(Np, 1)), R);
%             XX = 2 * rand(n, 1) - 1;
            f = 2 * liftFun(XX(:,k))' * A_bar' * Q_bar * Z;
            H = Z' * Q_bar * Z + R_bar;
            H = (H + H') / 2;
            if Constrained_MPC
    %             disp('Constrained On')
                u_seq = quadprog(2 * H, f, [], [], A_eq, b_eq, -control_bound * ones(Np, 1), control_bound * ones(Np, 1));
            end
        
            if Unconstrained_MPC
    %             disp('Unconstrained On')
                u_seq = -inv(2 * H) * f';
            end
            
            %% Weights Online Update
            U_seq = [U_seq u_seq(1:Nc)];
            X_seq = [X_seq A_bar * liftFun(XX(:,k)) + Z * u_seq];
        end
        E_D = norm(X_seq, 'fro');
        E_W = norm(U_seq, 'fro');
        diag_H = kron(eye(state_independent_num), 2 * H);
        gamma = N - 2 * r * (trace(        inv(diag_H)         )         );
        r = gamma / (2 * E_W);
        q = (nn - gamma) / (2 * E_D);
    
%         if gamma > 0
%             r = gamma / (2 * E_W);
%             q = (nn - gamma) / (2 * E_D);
%         end
        

        % Normalizing the number
        normalized_r = arrayfun(@(v) ceil(-log10(abs(v - round(v)))), r);
        max_decimal_digits_r = max(normalized_r);

        normalized_q = arrayfun(@(v) ceil(-log10(abs(v - round(v)))), q);
        max_decimal_digits_q = max(normalized_q);

        max_decimal_digits = max(max_decimal_digits_r, max_decimal_digits_q);

%         r = r * 10^max_decimal_digits
%         q = q * 10^max_decimal_digits

        r_set = [r_set r];
        q_set = [q_set q];
        q_r_relative_Set = [q_r_relative_Set q / r];
    end

end

% q = q * (1 / r)
% r = r * (1 / r)

X_Set = [x];
U_Set = [];
J_Set = [];
J = 0;
Online_Update = 0;

for i = 1 : 200
    if Online_Update
        for j = 1 : iterative_step
        
            Q = q * eye(n, n);
            R = r;
            Q_bar = kron(diag(ones(Np, 1)), Q);
            R_bar = kron(diag(ones(Np, 1)), R);
    
            Q_bar(Nc * n - n + 1 : Nc * n, Nc * n - n + 1 : Nc * n) = P;
            f = 2 * x' * A_bar' * Q_bar * Z;
            H = Z' * Q_bar * Z + R_bar;
            if Constrained_MPC
    %             disp('Constrained On')
                u_seq = quadprog(2 * H, f, [], [], A_eq, b_eq, -control_bound * ones(Np, 1), control_bound * ones(Np, 1));
            end
        
            if Unconstrained_MPC
    %             disp('Unconstrained On')
                u_seq = -inv(2 * H) * f';
            end
            %% Weights Online Update
            w = 0.2 * rand(n, 1) - 0.1;
            w = 0;
            x_seq = A_bar * x + Z * u_seq + w;
            E_D = x_seq' * x_seq;
            E_W = u_seq' * u_seq;
            gamma = N - 2 * r * (trace(2 * H))^(-1);
%             r = gamma / (2 * E_W)
%             q = (nn - gamma) / (2 * E_D)
            if E_D <= 1 || E_W <= 1
                r_set = [r_set r];
                q_set = [q_set q];
                continue;
            end
                
            if gamma > 0
                r = gamma / (2 * E_W)
                q = (nn - gamma) / (2 * E_D)
            end
            
            
            r_set = [r_set r];
            q_set = [q_set q];
        end
    end
    
    Q = q * eye(n, n);
    R = r;
    Q_bar = kron(diag(ones(Np, 1)), Q);
    R_bar = kron(diag(ones(Np, 1)), R);
    
    f = 2 * x' * A_bar' * Q_bar * Z;
    H = Z' * Q_bar * Z + R_bar;
    H = (H + H') / 2;
    if Constrained_MPC
%         disp('Constrained On')
        u_seq = quadprog(2 * H, f, [], [], A_eq, b_eq, -control_bound * ones(Np, 1), control_bound * ones(Np, 1));
    end
    
    if Unconstrained_MPC
%         disp('Unconstrained On')
        u_seq = -inv(2 * H) * f';
    end
    
    %% Continue dynamics
    w = 0.2 * rand(n, 1) - 0.1;
    w = 0;
%     x_next = A * x + B * u_seq(1) + w;

    x_next = f_ud(0, C * x, u_seq(1));
%     J = J + x' * Q * x + u_seq(1)' * R * u_seq(1)';
    x = liftFun(x_next);
    X_Set = [X_Set x];
    U_Set = [U_Set u_seq(1)];
%     J_Set = [J_Set J];

end
    

% figure 
% plot(r_set(1, :), 'LineStyle','-', 'LineWidth', 3.0) 
% legend('R')
% figure
% plot(q_set(1, :), 'LineStyle','-', 'LineWidth', 3.0) 
% legend('Q')
% figure
% plot(q_r_relative_Set(1,:), 'LineStyle','-', 'LineWidth', 3.0)
% legend('Q/R')
% 
% figure 
% plot(log(r_set(1, :)), 'LineStyle','-', 'LineWidth', 3.0) 
% legend('log(R)')
% figure
% plot(log(q_set(1, :)), 'LineStyle','-', 'LineWidth', 3.0) 
% legend('log(Q)')
% figure
% plot(log(q_r_relative_Set(1,:)), 'LineStyle','-', 'LineWidth', 3.0)
% legend('log(Q/R)')
% 
% 
% figure(1)
% plot(X_Set(1, :), 'LineStyle','-', 'LineWidth', 2.0) 
% legend('$x_1$ without Bayesian', '$x_1$ with Bayesian', 'Interpreter', 'latex')
% 
% figure(2)
% plot(X_Set(2, :), 'LineStyle','-', 'LineWidth', 2.0) 
% legend('$x_2$ without Bayesian', '$x_2$ with Bayesian', 'Interpreter', 'latex')
% 
% 
% 
% figure(3)
% plot(U_Set(1, :), 'LineStyle','-', 'LineWidth', 2.0) 
% hold on
% legend('$u$ without Bayesian', '$u$ with Bayesian', 'Interpreter', 'latex')
% 
% figure 
% plot(X_Set(1, :), 'LineStyle','-', 'LineWidth', 2.0) 
% hold on
% plot(X_Set(2, :), 'LineStyle','-', 'LineWidth', 2.0) 
% legend('$x_1$ with Bayesian', '$x_2$ with Bayesian', 'Interpreter', 'latex')
% 
% 
% 


% figure(4)
% plot(J_Set(1, :), 'LineStyle','-', 'LineWidth', 2.0)
% hold on 
% legend('$J$ without Bayesian', '$J$ with Bayesian', 'Interpreter', 'latex')
% 
% 


%% 画图
figure 
plot(r_set(1, :), 'LineStyle','-', 'LineWidth', 3.0) 
legend('$r$', 'Interpreter','latex')
% title('$r$', 'Interpreter','latex')
xlabel('Iterative steps', 'Interpreter','latex')
ylabel('$r$', 'Interpreter','latex')

figure
plot(q_set(1, :), 'LineStyle','-', 'LineWidth', 3.0) 
legend('$q$', 'Interpreter','latex')
% title('$q$', 'Interpreter','latex')
xlabel('Iterative steps', 'Interpreter','latex')
ylabel('$q$', 'Interpreter','latex')

figure(1)
plot(X_Set(1, :), 'LineStyle','-', 'LineWidth', 2.0) 
legend('$x_1$ with $q(0),r(0)$', '$x_1$ with $q^*,r^*$', 'Interpreter', 'latex')
xlabel('Iterative steps', 'Interpreter','latex')
ylabel('$x_1$', 'Interpreter','latex')


figure(2)
plot(X_Set(2, :), 'LineStyle','-', 'LineWidth', 2.0) 
legend('$x_2$ with $q(0),r(0)$', '$x_2$ with $q^*,r^*$', 'Interpreter', 'latex')
xlabel('Iterative steps', 'Interpreter','latex')
ylabel('$x_2$', 'Interpreter','latex')


figure(3)
plot(U_Set(1, :), 'LineStyle','-', 'LineWidth', 2.0) 
hold on
legend('$u$ with $q(0),r(0)$', '$u$ with $q^*,r^*$', 'Interpreter', 'latex')
xlabel('Iterative steps', 'Interpreter','latex')
ylabel('$u$', 'Interpreter','latex')

figure(12)
plot(q_set(1, :) ./ r_set(1,:), 'LineStyle','-', 'LineWidth', 3.0) 
legend('$q/r$', 'Interpreter','latex')
% title('$r$', 'Interpreter','latex')
xlabel('Iterative steps', 'Interpreter','latex')
ylabel('$q/r$', 'Interpreter','latex')


figure(10)
hold on
plot(X_Set(1, :), 'LineStyle','-', 'LineWidth', 3.0) 
hold on
plot(X_Set(2, :), 'LineStyle','-', 'LineWidth', 3.0) 
legend('$x_1$ with $q(0),r(0)$', '$x_2$ with $q(0),r(0)$', '$x_1$ with $q^*,r^*$', '$x_2$ with $q^*,r^*$', 'Interpreter', 'latex')
xlabel('Iterative steps', 'Interpreter','latex')
ylabel('$x$', 'Interpreter','latex')

