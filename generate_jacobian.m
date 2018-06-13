function A = generate_jacobian(f, x, h)
% Generates jacobian of vector function f using a step size of h
%   uses simple 3-point rule, O(h^2) accuracy
%   requires f be defined on a neighbourhood centered at x with a radius
%   larger than h

    % Do not permit a step greater than 10^-4 or less than 10^-11
    if (h > 1e-4 || h < 1e-11)
        % choice of sqrt(eps)*x obtained from:
        % https://en.wikipedia.org/wiki/Numerical_differentiation#
        %       Practical_considerations_using_floating_point_arithmetic
        h = sqrt(eps())*norm(x, inf);
    end
    
    m = size(f(x)); m = max(m(1),m(2));
    n = size(x); n = max(n(1),n(2));
    
    A = zeros(n, m);
    xp = x;
    xm = x;
    
    for i = 1:n
        xp(i) = xp(i) + h;
        xm(i) = xm(i) - h;
        
        A(:,i) = (f(xp)-f(xm))./(2*h);
        
        xp(i) = xp(i) - h;
        xm(i) = xm(i) + h;
    end
end