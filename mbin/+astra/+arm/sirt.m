function f = sirt(W,f_in,p,no_it)
%--------------------------------------------------------------------------
% f_new = sirt(W,f_in,p,no_it)
%
% Performs no_it SIRT iterations on the system Wf=p with initial guess f_in
%
% Input: W         The projection matrix; f_in      The initial guess; p
% The right-hand-side; no_it     Number of iterations to be performed.
%
% Output: f_new     The approximate solution for Wf=p after no_it
% iterations.
%--------------------------------------------------------------------------

global F_seg free

if isnumeric(W)
    %% matrix
    % Determine scaling matrices
    C = sum(W);
    R = sum(W,2);
    C(C==0) = 1;
    R(R==0) = 1;
    C = 1./C';
    R = 1./R;
    % C = diag(sparse(1./sum(W))); R = diag(sparse(1./sum(W,2)));
    
    % Remove appearances of inf
    C(isinf(C)) = 0;
    R(isinf(R)) = 0;
    
    % Start SIRT
    f = f_in;
    
    disp('SIRT');
    
    for i = 1:no_it
        r = p - W*f;
        
        %     f_new = f_old + (C*(W'*(R*r)));
        f = f + (C.*(W'*(R.*r)));
        
%         show(f);
%         imwrite(reshape(f,[256,256]), sprintf('/ufs/bleichro/Documents/Talks/Mathematisch_Congres/img/animproj/sirtrecon%d.png', i));
%         for j = 1:3
%             imwrite(reshape(f, [64,64]), sprintf('/ufs/bleichro/Documents/Talks/Mathematisch_Congres/img/animproj/simpel%d.png', (i-1)*3+j));
%         end
%         show(f);

        if false %~isempty(F_seg) && ~isempty(free)
            disp('Showing');
            % animate boundary update
            x = F_seg;
            x(free) = f;
            x(x>1) = 1;
            x(x<0) = 0;
            show(x, 'in', 'fit');
            for j = 1:6
                k = (i-1)*6+j;
                imwrite(x, sprintf('/ufs/bleichro/Documents/Talks/Mathematisch_Congres/img/animproj/update%d.png', k));
            end
        end
    end
    
elseif isa(W, 'opSpot')
    %% faster ASTRA implementation
    gpuEnabled = strcmpi(W.model_type, 'cuda');
    
    if gpuEnabled
        [~, f] = astra_create_reconstruction_cuda('SIRT_CUDA', W.proj_geom, W.vol_geom, reshape(p,W.proj_size), no_it, 'no', 0, 'no', 0, 'no', 0);
    else
        proj_id = astra_create_projector(W.model_type, W.proj_geom, W.vol_geom);
        [~, f] = astra_create_reconstruction('SIRT', proj_id, reshape(p,W.proj_size), no_it, 'no', 0, 'no', 0, 'no', 0);
    end
    f = f(:);
else
    error('Wrong matrix type!');
end
