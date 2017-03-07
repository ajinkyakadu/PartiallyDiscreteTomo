classdef opTV3D < opSpot
    % opTV, anisotropic TV operator
    properties (Access = private)
        funHandle
        hFilter = zeros(2,2,2);
        vFilter = zeros(2,2,2);
        zFilter = zeros(2,2,2);
        hFilterTranspose = zeros(3,3,3);
        vFilterTranspose = zeros(3,3,3);
        zFilterTranspose = zeros(3,3,3);
        imSize
    end
    
    methods
        
        % constructor
        function op = opTV3D(imSize)
            
            if numel(imSize) == 1
                imSize = [imSize, imSize, imSize];
            elseif numel(imSize) ~= 3
                error('Invalid image dimensions!');
            end
            
            N = prod(imSize);
            
            op = op@opSpot('opTV3D', 3*N, N);
            
            op.funHandle = @opTV_internal;
            
            op.imSize = imSize;
            
            op.hFilter(1,1:2,1) = [-1,1];
            op.vFilter(1:2,1,1) = [-1,1];
            op.zFilter(1,1,1:2) = [-1,1];
            
            op.hFilterTranspose(1,1:2,1) = [1,-1];
            op.vFilterTranspose(1:2,1,1) = [1,-1];
            op.zFilterTranspose(1,1,1:2) = [1,-1];
            
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
                y = zeros(3*op.n,1);
                
                % reshape x to image dimensions
                x = reshape(x, op.imSize);
                
                % horizontal filtering
                y(1:op.n) = imfilter(x, op.hFilter, 'symmetric');
                
                % vertical filtering
                y(op.n+1:2*op.n) = imfilter(x, op.vFilter, 'symmetric');
                
                % slice filtering
                y(2*op.n+1:end) = imfilter(x, op.zFilter, 'symmetric');
                
            % transposed product
            elseif mode == 2
                % TODO: do not waste memory!
                
                % reshape tv images
                tv_hor   = reshape(x(1:op.n), op.imSize);
                tv_ver   = reshape(x(op.n+1:2*op.n), op.imSize);
                tv_slice = reshape(x(2*op.n+1:end), op.imSize);
                
                y1 = imfilter(tv_hor, op.hFilterTranspose);
                y1(end,:,:) = tv_hor(end-1,:,:);
                
                y2 = imfilter(tv_ver, op.vFilterTranspose);
                y2(:,end,:) = tv_ver(:,end-1,:);
                
                y3 = imfilter(tv_slice, op.zFilterTranspose);
                y3(:,:,end) = tv_slice(:,:,end-1);
                
                y = y1(:) + y2(:) + y3(:);
            else
                % Should we really catch this?
                error('Multiplication error, mode should be 1 or 2.');
            end
            
        end
    end      
    
end
    