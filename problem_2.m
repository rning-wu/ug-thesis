g = 1;
n = 1;
F = @(t,v,x) -sin(x);
M = eye(n);
D = zeros(n);
K = eye(n);
t_start = 0;
t_end = 6*pi;
h = 0.1;
theta = 0.07;
rho_inf = 0.60;
x0 = 2;
v0 = 0;
[t_theta, x_theta, v_theta, a_theta, cond_theta, c_theta] = trap_bdf2_theta(F, M, D, K, x0, v0, h, t_start, t_end, theta);
[t_alpha, x_alpha, v_alpha, a_alpha, cond_alpha, c_alpha] = gen_alpha(F, M, D, K, x0, v0, h, t_start, t_end, rho_inf);

figure(20);
hold on;
plot(t_alpha, x_alpha);
plot(t_theta, x_theta);

legend('\alpha-method', '\theta-method', 'Location', 'Best')
fprintf('Number of nonlinear solves (bdf-theta): %i\n', c_theta);
fprintf('Number of nonlinear solves (gen-alpha): %i\n', c_alpha);