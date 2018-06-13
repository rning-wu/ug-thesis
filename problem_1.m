k = 10;
n = 2;
F = @(t,v,x) x-x;
M = eye(n);
D = zeros(n);
K = diag([1 k^2]);
t_start = 0;
t_end = 4980;
h = 0.1;
theta = 0.07;
rho_inf = 0.60;
x0 = [1 1]';
v0 = [0 0]';
[t_theta, x_theta, v_theta, a_theta, cond_theta, c_theta] = trap_bdf2_theta(F, M, D, K, x0, v0, h, t_start, t_end, theta);
[t_alpha, x_alpha, v_alpha, a_alpha, cond_alpha, c_alpha] = gen_alpha(F, M, D, K, x0, v0, h, t_start, t_end, rho_inf);

figure(10);
hold on;

plot(t_alpha, x_alpha(:,1));
plot(t_theta, x_theta(:,1));
xlabel('time')
ylabel('displacement')
legend('\alpha-method', '\theta-method', 'Location', 'Best')

figure(11);
hold on;

cut = 1500;

plot(t_alpha(1:cut), x_alpha(1:cut,2));
plot(t_theta(1:cut), x_theta(1:cut,2));
xlabel('time')
ylabel('displacement')
legend('\alpha-method', '\theta-method', 'Location', 'Best')

fprintf('Number of nonlinear solves (bdf-theta): %i\n', c_theta);
fprintf('Number of nonlinear solves (gen-alpha): %i\n', c_alpha);

figure(13);
hold on;
plot(1:length(cond_theta), cond_theta);
plot(1:length(cond_alpha), cond_alpha);
