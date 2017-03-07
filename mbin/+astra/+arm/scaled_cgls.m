function x_new = scaled_cgls(A,b,x_in,no_it)
%--------------------------------------------------------------------------
% x_new = scaled_cgls(A,b,x_in,no_it)
%               
% Performs no_it Scaled CGLS iterations on the system Wf=p with initial 
% guess f_in. The scaling corresponds to that of SIRT.
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

% Determine scaling matrices
C = diag(sparse(1./sum(A)));
R = diag(sparse(1./sum(A,2)));

% Remove appearances of inf
C(isinf(C)) = 0;
R(isinf(R)) = 0;

% Start scaled_CGLS
r_new = b - A*x_in;
z_tilde_new = A'*(R*r_new);
z_new = C*z_tilde_new;
p_new = z_new;
x_new = x_in;
inp_z_new = z_new'*z_tilde_new;

for i = 1:no_it

    x_old = x_new;
    r_old = r_new;    
    
    % Check if residual of the normal equations is or the preconditioned
    % residual is (almost) zero, if so, convergence is reached.
    if (norm(z_tilde_new,2) < eps || norm(z_new,2) < eps)
        display(['CGLS_reg converged in ' num2str(i) 'iterations.']);
        return
    end
    
    p_old = p_new;
    inp_z_old = inp_z_new;    
    
    w = A*p_old;
    alpha = (inp_z_old)/(w'*(R*w));
    x_new = x_old + alpha*p_old;
    r_new = r_old - alpha*w;
    z_tilde_new = A'*(R*r_new);
    z_new = C*z_tilde_new;
    inp_z_new = z_new'*z_tilde_new;
    beta = (inp_z_new)/(inp_z_old);
    p_new = z_new + beta*p_old;
    
end

end

