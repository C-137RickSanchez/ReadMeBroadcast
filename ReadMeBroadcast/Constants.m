classdef Constants
% This class defines the constants that are used for calculations
% and can be expanded with additional properties or methods for future applications
% Written by locateTempUserBash
% visit the user page @ github for further information
% or email using locateTempUserBash@yahoo.com
    properties (Constant)
        my=3.986005e14;             % Gravitional constant for WGS84
        OmegaDotE=7.2921151467e-5;  % Earth's rotation rate for WGS84             
        C=299792458;                % The speed of light in vacuum            
         
    end
    methods(Static)
     function result = getConstant(prop)
                    result = obj.(prop);
     end
    end
end