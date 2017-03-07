classdef algorithm_handle < handle
    %ASTRA.ALGORITHM_HANDLE Handle class around an astra_mex_algorithm id
    %   Automatically deletes the algorithm when deleted.
    
    properties
        id
    end
    
    methods
        function obj = algorithm_handle(alg_id)
            obj.id = alg_id;
        end
        function delete(obj)
            astra_mex_algorithm('delete', obj.id);
        end
    end
    
end

