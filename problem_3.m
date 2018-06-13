h = 0.05;
its = 3;
epsi = 1e-2; % to be paramatrized

r = 1;

F = @(t, v, x) [
    -epsi^(-2)*(x(1) - 1) + x(1)^(-3)*(v(2)^2) + cos(x(2));
    -x(1)^(-1)*sin(x(2))
    ];
M = eye(2);
D = zeros(2);
K = zeros(2);

x0 = [1+epsi, pi/4]';
v0 = [1/sqrt(2), -1/sqrt(2)]';

y0 = [r+epsi; 1/sqrt(2); pi/4; -1/sqrt(2)];
tr = [0 16*pi];

figure(31);
theta = 0.07;
[t_theta, x_theta, v_theta, ~, cond_theta, c_theta] = trap_bdf2_theta(F, M, D, K, x0, v0, h, tr(1), tr(2), theta);
plot(t_theta, x_theta);
legend('High frequency', 'Low frequency', 'Location', 'Best')
xlabel('time')
ylabel('displacement')

figure(32);
rho_inf = 0.6;
[t_alpha, x_alpha, v_alpha, ~, cond_alpha, c_alpha] = gen_alpha(F, M, D, K, x0, v0, h, tr(1), tr(2), rho_inf);
plot(t_alpha, x_alpha);
legend('High frequency', 'Low frequency', 'Location', 'Best')
xlabel('time')
ylabel('displacement')

figure(33);
hold on;
plot(1:length(cond_alpha), cond_alpha, 'b.-');
plot(1:length(cond_theta), cond_theta, 'r.-');

legend('\alpha-method', '\theta-method', 'Location', 'Best')
xlabel('time')
ylabel('\kappa')


fprintf('Number of nonlinear solves (bdf-theta): %i\n', c_theta);
fprintf('Number of nonlinear solves (gen-alpha): %i\n', c_alpha);