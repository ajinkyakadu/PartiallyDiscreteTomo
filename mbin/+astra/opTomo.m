% -----------------------------------------------------------------------
% Copyright 2014 Centrum Wiskunde & Informatica, Amsterdam
%
% Author:  Folkert Bleichrodt
% Contact: F.Bleichrodt@cwi.nl
% 
% This file is the ASTRA-Spot operator "opTomo" to be used with the
% All Scale Tomographic Reconstruction Antwerp Toolbox ("ASTRA Toolbox").
%
% The ASTRA-Spot operator is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The ASTRA-Spot operator is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the ASTRA-Spot operator. If not, see <http://www.gnu.org/licenses/>.
%
%-----------------------------------------------------------------------
classdef opTomo < opSpot
    %OpTomo  Wrapper for Astra tomography projector
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = private )
        % Multiplication function
        funHandle
        % ASTRA identifiers
        sino_id
        vol_id
        fp_alg_id
        bp_alg_id
        % ASTRA IDs handle
        astra_handle
        % mask
        mask = [];
    end % Properties

    properties ( SetAccess = private, GetAccess = public )
        proj_size
        vol_size
    end % properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % opTomo. Constructor.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opTomo(type, proj_geom, vol_geom, gpu_index)
            
            if nargin < 4 || isempty(gpu_index), gpu_index = 0; end
            
            proj_size = astra_geom_size(proj_geom);
            vol_size  = astra_geom_size(vol_geom);

            % construct operator
            op = op@opSpot('opTomo', prod(proj_size), prod(vol_size));
            
            % determine the dimension
            % TODO: Check consistency between proj_geom and vol_geom
            is2D = ~isfield(vol_geom, 'GridSliceCount');
            gpuEnabled = strcmpi(type, 'cuda');
            
            if is2D
                % create a projector
                proj_id = astra_create_projector(type, proj_geom, vol_geom);
                
                % create a function handle
                op.funHandle = @opTomo_intrnl2D;
                
                % Initialize ASTRA data objects.
                % projection data
                sino_id = astra_mex_data2d('create', '-sino', proj_geom, 0);

                % image data
                vol_id  = astra_mex_data2d('create', '-vol', vol_geom, 0);
                
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
                % volume data
                vol_id  = astra_mex_data3d('create', '-vol', vol_geom, 0);
                
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
            
            fp_alg_id = astra_mex_algorithm('create', cfg_fp);
            bp_alg_id = astra_mex_algorithm('create', cfg_bp);
            
            % Create handle to ASTRA objects, so they will be deleted if
            % opTomo is deleted.
            op.astra_handle = astra.opTomo_helper_handle([sino_id, vol_id, ...
                proj_id, fp_alg_id, bp_alg_id]);
            
            % pass object properties
            op.proj_size   = proj_size;
            op.vol_size    = vol_size;          
%             op.model_type  = type;
%             op.proj_geom   = proj_geom;
%             op.vol_geom    = vol_geom;
            op.fp_alg_id   = fp_alg_id;
            op.bp_alg_id   = bp_alg_id;
            op.sino_id     = sino_id;
            op.vol_id      = vol_id;
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
            
            % ASTRA cannot handle sparse vectors
            if issparse(x)
                x = full(x);
            end
            
            if mode == 1              
                % X is passed as a vector. Reshape it into an image.             
                x = reshape(x, op.vol_size);
                
                % Matlab data copied to ASTRA data
                astra_mex_data2d('store', op.vol_id, x);
                
                % forward projection
                astra_mex_algorithm('iterate', op.fp_alg_id);
                
                % retrieve Matlab array. Reshape it back into a vector
                y = reshape(astra_mex_data2d('get', op.sino_id),[],1);
            else
                % X is passed as a vector. Reshape it into a sinogram.
                x = reshape(x, op.proj_size);
                
                % Matlab data copied to ASTRA data
                astra_mex_data2d('store', op.sino_id, x);
                
                % backprojection
                astra_mex_algorithm('iterate', op.bp_alg_id);
                
                % retrieve Matlab array. Reshape it back into a vector
                y = reshape(astra_mex_data2d('get', op.vol_id),[],1);
            end
        end


        % 3D projection code
        function y = opTomo_intrnl3D(op,x,mode)
            
            % ASTRA cannot handle sparse vectors
            if issparse(x)
                x = full(x);
            end
            
            if mode == 1
                % X is passed as a vector. Reshape it into an image.
                x = reshape(x, op.vol_size);

                % Matlab data copied to ASTRA data
                astra_mex_data3d('store', op.vol_id, x);
                
                % forward projection
                astra_mex_algorithm('iterate', op.fp_alg_id);
                
                % retrieve Matlab array. Reshape it back into a vector
                y = reshape(astra_mex_data3d('get', op.sino_id),[],1);
            else
                % X is passed as a vector. Reshape it into projection data.
                x = reshape(x, op.proj_size);
                
                % Matlab data copied to ASTRA data
                astra_mex_data3d('store', op.sino_id, x);

                % forward projection
                astra_mex_algorithm('iterate', op.bp_alg_id);
                
                % retrieve Matlab array. Reshape it back into a vector
                y = reshape(astra_mex_data3d('get', op.vol_id),[],1);
            end
        end
    end % Methods
 
end % classdef
