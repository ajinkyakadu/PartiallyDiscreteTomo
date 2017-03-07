classdef optomo_helper_handle < handle
    %ASTRA.OPTOMO_HELPER_HANDLE Handle class around an astra identifier
    %   Automatically deletes the data when deleted. Multiple id's can be
    %   passed as an array as input to the constructor.
    
    properties
        id
    end
    
    methods
        function obj = astra_handle(id)
            obj.id = id;
        end
        function delete(obj)
            for i = 1:numel(obj.id)
                astra_mex_data2d('delete', obj.id(i));
                astra_mex_data3d('delete', obj.id(i));
                astra_mex_algorithm('delete', obj.id(i));
                astra_mex_matrix('delete', obj.id(i));
                astra_mex_projector('delete', obj.id(i));
                astra_mex_projector3d('delete', obj.id(i))
            end
        end
    end
    
end

