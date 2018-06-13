% Takes one step 

function [v1, iterations, condition] = be_trap_step(F, M, D, K, x0, v0, a0, h, t, theta)
    theta_original = theta;
    theta = 0.50 + theta/50;
    
    f = @(v) v - (v0 + h*((1-theta)*a0 + theta*(M\F(t, v0, x0) - D*v - K*(x0 + h*((1-theta)*v0 + theta*v)))));
        
    [v1, iterations, success] = nlfzero(f, v0, 3);
    
    if (~success && h > 1e-5)
        [vm, it1, ~] = be_trap_step(F, M, D, K, x0, v0, a0, h/2, t, theta);        
        xm = x0 + h*((1-theta)*v0 + theta*vm);
        am = M\(F(t+h/2, vm, xm) - D*vm - K*xm);
        
        theta = theta_original;
        
        [v1, it2, ~] = trap_bdf2_theta_step(F, M, D, K, x0, xm, v0, vm, am, h/2, t+h/2, theta);
        iterations = iterations + it1 + it2;
    end
    
    condition = cond(generate_jacobian(f, a0, 1e-6));
end