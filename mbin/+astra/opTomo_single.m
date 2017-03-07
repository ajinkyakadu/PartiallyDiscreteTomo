classdef opTomo_single < opSpot
    %OpTomo  Wrapper for Astra tomography projector

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = private )
        funHandle       % Multiplication function
        proj_handle
        % ASTRA data identifiers
        sino_id
        vol_id
        sino_handle
        vol_handle
        % ASTRA algorithm identifiers
        fp_alg_id
        bp_alg_id
        % mask
        mask = [];
    end % Properties

    properties ( SetAccess = private, GetAccess = public )
        proj_size
        vol_size
        proj_geom
        vol_geom
        model_type
    end % properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % opTomo. Constructor.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opTomo_single(type, proj_geom, vol_geom, gpu_index)
            
            if nargin < 4 || isempty(gpu_index), gpu_index = 0; end
            
            proj_size = astra_geom_size(proj_geom);
            vol_size  = astra_geom_size(vol_geom);
            proj_size = astra_geom_size(proj_geom);
            vol_size  = astra_geom_size(vol_geom);

            % construct operator
            op = op@opSpot('opTomo_single', prod(proj_size), prod(vol_size));

            
            % determine the dimension
            % TODO: Check consistency between proj_geom and vol_geom
            is2D = ~isfield(vol_geom, 'GridSliceCount');
            gpuEnabled = strcmpi(type, 'cuda');
            
            if is2D
                % create a projector
                proj_id = astra_create_projector(type, proj_geom, vol_geom);
                % handle, for freeing projector if object is deleted
                proj = astra.projector_handle(proj_id);
                
                % create a function handle
                op.funHandle = @opTomo_intrnl2D;
                
                % Initialize ASTRA data objects.
                % projection data
                sino_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
                % handle, for freeing data if object is deleted
                sino_handle = astra.data2d_handle(sino_id);
                % image data
                vol_id  = astra_mex_data2d('create', '-vol', vol_geom, 0);
                % handle, for freeing data if object is deleted
                vol_handle = astra.data2d_handle(vol_id);
                
                % forward and back projection algorithms
                if gpuEnabled
                    fp_alg = 'FP_CUDA';
                    bp_alg = 'BP_CUDA';
                    proj_id = [];
                else
                    fp_alg = 'FP';
                    bp_alg = 'BP';
                    proj_id = astra_create_projector(type, proj_geom, vol_geom);
                end
            else
                % only gpu/cuda code for 3d
                if ~gpuEnabled
                    error(['Only type ' 39 'cuda' 39 ' is supported ' ...
                           'for 3D geometries.'])
                end
                
                % projector is builtin for cuda code
                proj = [];
                proj_id = [];
                
                % create a function handle
                op.funHandle = @opTomo_intrnl3D;
                
                % Initialize ASTRA data objects.
                % projection data
                sino_id = astra_mex_data3d('create', '-sino', proj_geom, 0);
                % handle, for freeing data if object is deleted
                sino_handle = astra.data3d_handle(sino_id);
                % volume data
                vol_id  = astra_mex_data3d('create', '-vol', vol_geom, 0);
                % handle, for freeing data if object is deleted
                vol_handle = astra.data3d_handle(vol_id);
                
                % forward and back projection algorithms
                fp_alg = 'FP3D_CUDA';
                bp_alg = 'BP3D_CUDA';
            end
 
            % Configuration for astra fp algorithm
            cfg_fp = astra_struct(fp_alg);
            cfg_fp.ProjectorId      = proj_id;
            cfg_fp.ProjectionDataId = sino_id;
            cfg_fp.VolumeDataId     = vol_id;
            if gpuEnabled
                cfg_fp.option.GPUindex  = gpu_index;
            end

            cfg_bp = astra_struct(bp_alg);
            cfg_bp.ProjectionDataId     = sino_id;
            cfg_bp.ProjectorId          = proj_id;
            cfg_bp.ReconstructionDataId = vol_id;
            if gpuEnabled
                cfg_bp.option.GPUindex      = gpu_index;
            end
            
            % pass object properties
            op.proj_size   = proj_size;
            op.vol_size    = vol_size;          
            op.model_type  = type;
            op.proj_handle = proj;
            op.proj_geom   = proj_geom;
            op.vol_geom    = vol_geom;
            op.fp_alg_id   = astra_mex_algorithm('create', cfg_fp);
            op.bp_alg_id   = astra_mex_algorithm('create', cfg_bp);
            op.sino_id     = sino_id;
            op.vol_id      = vol_id;
            op.sino_handle = sino_handle;
            op.vol_handle  = vol_handle;
            op.cflag       = false;
            op.sweepflag   = false;

        end %

    end % methods - public

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )

        % Multiplication
        function y = multiply(op,x,mode)
            y = op.funHandle(op, x, mode);
        end % Multiply

    end % Methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - private
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = private )

        % 2D projection code
        function y = opTomo_intrnl2D(op,x,mode)
            if mode == 1
                % X is passed as a vector. Reshape it into an image.
                x = reshape(x, op.vol_size);
                
                % Matlab data copied to ASTRA data
                astra_mex_data2d('store', op.vol_id, x);
                
                % forward projection
                astra_mex_algorithm('iterate', op.fp_alg_id);
                
                % retrieve Matlab array. Reshape it back into a vector
                y = reshape(astra_mex_data2d('get_single', op.sino_id),[],1);
            else
                % X is passed as a vector. Reshape it into a sinogram.
                x = reshape(x, op.proj_size);
                
                % Matlab data copied to ASTRA data
                astra_mex_data2d('store', op.sino_id, x);
                
                % backprojection
                astra_mex_algorithm('iterate', op.bp_alg_id);
                
                % retrieve Matlab array. Reshape it back into a vector
                y = reshape(astra_mex_data2d('get_single', op.vol_id),[],1);
            end
        end


        % 3D projection code
        function y = opTomo_intrnl3D(op,x,mode)
            if mode == 1
                % X is passed as a vector. Reshape it into an image.
                x = reshape(x, op.vol_size);

                % Matlab data copied to ASTRA data
                astra_mex_data3d('store', op.vol_id, x);
                
                % forward projection
                astra_mex_algorithm('iterate', op.fp_alg_id);
                
                % retrieve Matlab array. Reshape it back into a vector
                y = reshape(astra_mex_data3d('get_single', op.sino_id),[],1);
            else
                % X is passed as a vector. Reshape it into projection data.
                x = reshape(x, op.proj_size);
                
                % Matlab data copied to ASTRA data
                astra_mex_data3d('store', op.sino_id, x);

                % forward projection
                astra_mex_algorithm('iterate', op.bp_alg_id);
                
                % retrieve Matlab array. Reshape it back into a vector
                y = reshape(astra_mex_data3d('get_single', op.vol_id),[],1);
            end
        end
    end % Methods
 
end % classdef
