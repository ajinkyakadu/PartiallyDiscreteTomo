function x_new = cgne(A,b,x_in,no_it)
%--------------------------------------------------------------------------
% x_new = cgne(A,b,x_in,no_it)
%               
% Performs no_it CGNE iterations on the system Wf=p with initial guess
% f_in. This method is known to be unstable when b contains noise.
%
% Input:
% A         The projection matrix;
% b         The right-hand-side;
% x_in      The initial guess;
% no_it     Number of iterations to be performed.
%
% Output:
% x_new     The approximate solution for Wf=p after no_it iterations.
%--------------------------------------------------------------------------

% Start CGNE
r_new = b - A*x_in;
p_new = A'*r_new;
x_new = x_in;

for i = 1:no_it

    x_old = x_new;
    r_old = r_new;
    
    % Check if residual is (almost) zero, if so, convergence is reached.
    if (norm(r_new,2) < eps)
        display(['CGNE converged in ' num2str(i) 'iterations.']);
        return
    end
    p_old = p_new;
    
    alpha = (r_old'*r_old)/(p_old'*p_old);
    x_new = x_old + alpha*p_old;
    r_new = r_old - alpha*A*p_old;
    beta = (r_new'*r_new)/(r_old'*r_old);
    p_new = A'*r_new + beta*p_old;
    
end

end

