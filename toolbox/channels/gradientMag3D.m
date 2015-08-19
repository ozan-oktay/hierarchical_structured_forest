function [M,O,varargout] = gradientMag3D( I, normRad, normConst )

residual = 1e-12;
if(nargin<1 || isempty(I)), M=single([]); O=M; return; end
if(nargin<2 || isempty(normRad)), normRad=0; end
if(nargin<3 || isempty(normConst)), normConst=.005; end

if (ndims(I)==2)  %#ok<*ISMAT>
    
    [Gx,Gy]  = gradient(I,1);
    Gy       = Gy + residual;
    M        = sqrt(Gx.^2+Gy.^2); 
    O        = atan2(Gy,Gx);
    O(O<0)   = O(O<0) + pi;

elseif (ndims(I)==3)

    [Gx,Gy,Gz]  = gradient(I,1);
    Gy          = Gy + residual;
    M           = sqrt(Gx.^2+Gy.^2+Gz.^2); 
    O           = atan2(Gy,Gx);
    O(O<0)      = O(O<0) + pi;  
        
else
    error('gradientMag3D::ndims error')
end

S = convTri3D( M, normRad );
M = M./(S + normConst);

if (nargout>=3 && ndims(I)==2)
    varargout{1} = cat(3,Gx,Gy);
elseif (nargout>=3 && ndims(I)==3)  
    varargout{1} = cat(4,Gx,Gy,Gz);
end

end

