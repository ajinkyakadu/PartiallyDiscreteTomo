classdef data3d_handle < handle
    %ASTRA.DATA3D_HANDLE Handle class around an astra_mex_data3d id
    %   Automatically deletes the data when deleted.
    
    properties
        id
    end
    
    methods
        function obj = data3d_handle(data_id)
            obj.id = data_id;
        end
        function delete(obj)
            astra_mex_data3d('delete', obj.id);
        end
    end
    
end

