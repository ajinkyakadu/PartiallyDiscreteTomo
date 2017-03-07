classdef data2d_handle < handle
    %ASTRA.DATA2D_HANDLE Handle class around an astra_mex_data2d id
    %   Automatically deletes the data when deleted.
    
    properties
        id
    end
    
    methods
        function obj = data2d_handle(data_id)
            obj.id = data_id;
        end
        function delete(obj)
            astra_mex_data2d('delete', obj.id);
        end
    end
    
end

