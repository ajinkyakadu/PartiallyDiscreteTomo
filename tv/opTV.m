classdef opTV < opSpot
    % opTV, anisotropic TV operator
    properties (Access = private)
        funHandle
        filter = [-1,1];
        filterTranspose = [1,-1,0]; 
        width
    end
    
    methods
        
        % constructor
        function op = opTV(n)
            % create operator
            
            N = n*n;
            
            op = op@opSpot('opTV', 2*N, N);
            
            op.funHandle = @opTV_internal;
            
            op.width = n;
            
            op.cflag = false;
            op.sweepflag = false;
        end
    end % Public methods
    
    methods(Access = protected)
        
        % Multiplication
        function y = multiply(op, x, mode)
            y = op.funHandle(op, x, mode);
        end
        
    end % Protected methods
    
    methods(Access = private)
        
        function y = opTV_internal(op, x, mode)
            
            % forward product
            if mode == 1
                
                % allocate output
                y = zeros(2*op.n,1);
                
                % reshape x to image dimensions
                x = reshape(x, [op.width,op.width]);
                
                % apply filtering (vertical)
                y(1:op.n) = imfilter(x, op.filter','symmetric');
                
                % apply filtering (horizontal)
                y(op.n+1:end) = imfilter(x, op.filter,'symmetric');
                
            % transposed product
            elseif mode == 2
                % TODO: do not waste memory!
                
                % reshape horizontal TV
                tv_hor = reshape(x(1:op.n), [op.width, op.width]);
                
                % reshape vertical TV
                tv_ver = reshape(x(op.n+1:end), [op.width, op.width]);
                
                y1 = imfilter(tv_hor, op.filterTranspose');
                y1(end,:) = tv_hor(end-1,:);
                
                y2 = imfilter(tv_ver, op.filterTranspose);
                y2(:,end) = tv_ver(:,end-1);
                
                y = y1(:) + y2(:);
            else
                % Should we really catch this?
                error('Multiplication error, mode should be 1 or 2.');
            end
            
        end
    end      
    
end
    