function [v2, iterations, condition] = trap_bdf2_theta_step(F, M, D, K, x0, x1, v0, v1, a1, h, t, theta)
    f = @(v) v - (v1 + (1-theta)*(h/2)*(((M\(F(t, v, x1 + (1-theta)*(h/2)*(v + v1) + (theta/3)*(x1-x0+2*h*v)) - D*v - K*(x1 + (1-theta)*(h/2)*(v + v1) + (theta/3)*(x1-x0+2*h*v)))) + a1)) + (theta/3)*(v1-v0+2*h*((M\(F(t, v, x1 + (1-theta)*(h/2)*(v + v1) + (theta/3)*(x1-x0+2*h*v)) - D*v - K*(x1 + (1-theta)*(h/2)*(v + v1) + (theta/3)*(x1-x0+2*h*v)))))));
    [v2, iterations, success] = nlfzero(f, v1, 3);
    
    if (~success && h > 1e-5)
        [ah, it, ~]= be_trap_step(F, M, D, K, x1, v1, a1, h/2, t, theta);
    
        vh = v1 + h/2*(theta*ah+(1-theta)*a1);
        xh = x1 + h/2*(theta*vh+(1-theta)*v1);
            
        [v2, it2, ~] = trap_bdf2_theta_step(F, M, D, K, x1, xh, v1, vh, ah, h, t, theta);
        iterations = iterations + it + it2;
    end
    
    condition = cond(generate_jacobian(f, v1, 1e-6));
end

