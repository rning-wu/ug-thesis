function [acc, iterations, condition] = gen_alpha_step(F, M, D, K, t, x0, v0, a0, h, rho_inf, steps)    
    alpha_m = (2*rho_inf-1)/(rho_inf+1);
    alpha_f = (rho_inf)/(rho_inf+1);
    gamma=1/2-alpha_m+alpha_f;
    beta=(1/4)+(alpha_f-alpha_m)/2;
    
    t1 = t + h;
    t_wt = (1-alpha_f)*t1 + alpha_f*t;
    
    iterations = 0;
    
    f = @(w) (M*((1-alpha_m)*w + alpha_m*a0) + ...
        D*((1-alpha_f)*(v0 + h*((1-gamma)*a0 + gamma*w)) + alpha_f*v0) + ...
        K*((1-alpha_f)*(x0 + h*v0 + h^2*((1/2 - beta)*a0 + beta*w)) + alpha_f*x0) - ...
        F(t_wt, ...
        (1-alpha_f)*(v0 + h*((1-gamma)*a0 + gamma*w)) + alpha_f*v0, ...
        (1-alpha_f)*(x0 + h*v0 + h^2*((1/2 - beta)*a0 + beta*w)) + alpha_f*x0));
    
    [acc, iter, success] = nlfzero(f, a0, steps);
    
    if (~success && h > 1e-3)
        [a1, it, ~]= gen_alpha_step(F, M, D, K, t, x0, v0, a0, h/2, rho_inf, steps);
        
        v1 = v0 + h*((1-gamma)*a0 + gamma*a1);            
        x1 = x0 + h*v0 + h^2*((1/2-beta)*a0+beta*a1); 
        
        [acc, it2, ~] = gen_alpha_step(F, M, D, K, t+h/2, x1, v1, a1, h/2, rho_inf, steps);
        iterations = iterations + it + it2;
    end
    iterations = iterations + iter;
    condition = cond(generate_jacobian(f, a0, 1e-6));
end