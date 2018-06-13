function [tS, xS, vS, aS, conditions, count] = trap_bdf2_theta(F, M, D, K, x0, v0, h, t_start, t_end, theta)
    N = int64(floor((t_end-t_start)/h))+1;

    tS = zeros(N, 1);
    xS = zeros(N, numel(x0));
    vS = zeros(N, numel(x0));
    aS = zeros(N, numel(x0));
    
    a0 = M\(F(t_start, v0, x0)-D*v0-K*x0);

    conditions = zeros(N-1, 1);
    
    [a1, count, conditions(1)] = be_trap_step(F, M, D, K, x0, v0, a0, h, t_start, 0.5);
    
    tS(1) = t_start;
    tS(2) = tS(1)+h;

    aS(1,:) = a0;
    aS(2,:) = a1;
    
    vS(1,:) = v0;
    vS(2,:) = v0 + h*(0.5*aS(2,:)'+0.5*aS(1,:)');
    
    xS(1,:) = x0;
    xS(2,:) = x0 + h*(0.5*vS(2,:)'+0.5*vS(1,:)');
    
    for i = 3:N
        tS(i) = tS(i-1) + h;
        [vS(i,:), c, cond] = trap_bdf2_theta_step(F, M, D, K, xS(i-2,:)', xS(i-1,:)', vS(i-2,:)', vS(i-1,:)', aS(i-1,:)', h, tS(i), theta);
        xS(i,:) = xS(i-1,:) + (1-theta)*h/2*(vS(i-1,:) + vS(i,:)) + theta/3*(xS(i-1,:) - xS(i-2,:) + 2*h*vS(i,:));
        aS(i,:) = M\(F(tS(i), vS(i,:)', xS(i,:)') - D*vS(i,:)' - K*xS(i,:)');
                
        conditions(i-1) = cond;
        count = count+c;
    end
end

