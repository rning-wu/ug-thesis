function [tS, xS, vS, aS, conditions, count] = gen_alpha(F, M, D, K, x0, v0, h, t_start, t_end, rho_inf)
% gen_alpha Generalized Alpha method of Chung & Hulbert
%   Generalized Alpha method of Chung & Hulbert
%   Parameters: F is a vector function of applied forces  (mx1)
%               M is the mass matrix                      (mxm)
%               D is the damping matrix                   (mxm)
%               K is the spring constant matrix           (mxm)
%               x0, v0 are the initial conditions         (mx1)
%               h is the step size                        (real number)
%               t_start is the start time                 (real number)
%               t_end is the end time                     (real number)
%               rho_inf, alpha_f, alpha_m, gamma, beta are as specified in
%               Chung & Hulbert, 1993                     (real number)

    alpha_m = (2*rho_inf-1)/(rho_inf+1);
    alpha_f = (rho_inf)/(rho_inf+1);
    gamma=1/2-alpha_m+alpha_f;
    beta=(1/4)+(alpha_f-alpha_m)/2;
    
    N = int64(floor((t_end-t_start)/h))+1;
    
    a0 = M\(F(t_start, x0, v0)-D*v0-K*x0);

    tS = zeros(N, 1);
    tS(1) = t_start;    
    xS = zeros(N, numel(x0));
    vS = zeros(N, numel(x0));
    aS = zeros(N, numel(x0));
    
    aS(1,:) = a0;
    vS(1,:) = v0;    
    xS(1,:) = x0;
    
    max_its = 3;
    conditions = zeros(N-1, 1);
    count = 0;
        
    for i = 1:N-1
        tS(i+1) = tS(i) + h;
        [aS(i+1,:),c,conditions(i)] = gen_alpha_step(F, M, D, K, tS(i), xS(i,:)', vS(i,:)', aS(i,:)', h, rho_inf, max_its);
        count = count + c;

        vS(i+1,:) = vS(i,:) + h*((1-gamma)*aS(i,:) + gamma*aS(i+1,:));            
        xS(i+1,:) = xS(i,:) + h*vS(i,:) + h^2*((1/2-beta)*aS(i,:)+beta*aS(i+1,:)); 
    end
end