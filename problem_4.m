% simulations for coupled equations
h = 0.05;
its = 3;
epsi = 1e-4; % to be paramatrized

r = 1;

F = @(t, v, x) [
    -epsi^(-2)*(sqrt(x(1)^2+x(2)^2)-r)*(x(1)/sqrt(x(1)^2+x(2)^2))-cos((x(2)/sqrt(x(1)^2+x(2)^2)).^2)*(x(1)/sqrt(x(1)^2+x(2)^2)); 
    -epsi^(-2)*(sqrt(x(1)^2+x(2)^2)-r)*(x(2)/sqrt(x(1)^2+x(2)^2))-cos((-x(1)/sqrt(x(1)^2+x(2)^2)).^2)*(x(2)/sqrt(x(1)^2+x(2)^2));
];

M = eye(2);
D = zeros(2);
K = zeros(2);
x0 = [1/sqrt(2)+epsi; 1/sqrt(2)-epsi];
v0 = [1+epsi; -epsi];

y0 = [1/sqrt(2)+epsi; 1/sqrt(2)-epsi; 1+epsi; -epsi];
tr = [0 16*pi];

% bdf-integration
figure(41);
hold on;
theta = 0.07;
[t_theta, x_theta, v_theta, a_theta, cond_theta, c_theta] = trap_bdf2_theta(F, M, D, K, x0, v0, h, tr(1), tr(2), theta);
plot(t_theta, sqrt(x_theta(:,1).^2 + x_theta(:,2).^2));
plot(t_theta, x_theta(:,1));
plot(t_theta, x_theta(:,2));
xlabel('time')
ylabel('displacement')
legend('radius', 'x_1', 'x_2', 'Location', 'Best')

% gen-alpha integration
figure(42);
hold on;
rho_inf = 0.6;
[t_alpha, x_alpha, v_alpha, a_alpha, cond_alpha, c_alpha] = gen_alpha(F, M, D, K, x0, v0, h, tr(1), tr(2), rho_inf);
plot(t_alpha, sqrt(x_alpha(:,1).^2 + x_alpha(:,2).^2));
plot(t_alpha, x_alpha(:,1));
plot(t_alpha, x_alpha(:,2));
xlabel('time')
ylabel('displacement')
legend('radius', 'x_1', 'x_2', 'Location', 'Best')

fprintf('Number of nonlinear solves (bdf-theta): %i\n', c_theta);
fprintf('Number of nonlinear solves (gen-alpha): %i\n', c_alpha);

figure(43);
hold on;
plot(1:length(cond_theta), cond_theta, 'r.-');
plot(1:length(cond_alpha), cond_alpha, 'b.-');
legend('\kappa for \theta method', '\kappa for gen-\alpha method', 'Location', 'Best')
xlabel('time')
ylabel('\kappa')
