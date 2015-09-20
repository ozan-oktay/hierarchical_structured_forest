function M = gradientMag3D( I, normRad, normConst )

if(nargin<1 || isempty(I)), M=single([]); return; end
if(nargin<2 || isempty(normRad)), normRad=0; end
if(nargin<3 || isempty(normConst)), normConst=.005; end

if (ndims(I)==2)  %#ok<*ISMAT>
    
    [Gx,Gy]  = gradient(I,1);
    M        = sqrt(Gx.^2+Gy.^2); 

elseif (ndims(I)==3)

    [Gx,Gy,Gz]  = gradient(I,1);
    M           = sqrt(Gx.^2+Gy.^2+Gz.^2); 
        
else
    error('gradientMag3D::ndims error')
end

S = convTri3D( M, normRad );
M = M./(S + normConst);

end

