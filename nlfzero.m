function [r, i, success] = nlfzero(f, x, iter)
% Attempts to find a zero of the vector-valued function f with x as an
% initial estimate. If the residual doesn't decrease after iter iterations,
% consider the attempt to have failed and the function which called this is
% responsible for addressing this issue. In the context of ODE problems,
% this is halving the time-step. 

% Output variables: 
% r: the vector which makes f(r) = 0 (or, reasonably close within 1e-6)
% i: the number of iterations it took to get from x to r
% success: logical value; is true if and only if r was successfully found
%       in less than 40 iterations, and if the residual didn't increase
%       after iter iterations. 

    success = true;
    r = x;
    i = 0;
    
    % 40 iterations has been set to the maximum for now. 
    while (i < 40)
        % If the 2-norm of the resudual is less than 1e-6, then consider
        % the iteration to have succeeded. 
        L = 2;
        % TODO: consider possibly replacing this by the max-norm and
        % relaxing the tolerance level. 
        if (norm(f(r), L) < 1e-6)
            success = true;
            break
        end
        
        % If the step has used more than iter iterations, and the residual
        % has not decreased (to within machine tolerance of 1e-8), consider
        % the step to have failed and substep. 
        if (i >= iter && norm(f(r), L) >= norm(f(x), L) - 1e-8)
            success = false;
            break;
        end
        J = generate_jacobian(f, r, 1e-6);
        if (rank(J,1e-12) < size(J,1))
            success = false;
            break;
        end
        r = r + J\(-f(r));
        i = i + 1;
    end
    
    % If it took us 40 iterations to get here, consider the iteration to 
    % have failed. The process which called this is responsible for
    % substepping. 
    if (success) 
        success = (i < 40);
    end 
end