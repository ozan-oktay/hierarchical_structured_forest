function mat_rs = imResample3D( varargin )

%RESIZE     Resize a matrix.

% DESCRIPTION:
%       Resize a matrix to a given size using interp2 (2D) or interp3
%       (3D).
%       Use interpolation to redivide the [0,1] interval into Nx, Ny, Nz 
%       voxels, where 0 is the center of first voxel, and 1 is the center 
%       of the last one.
%
% USAGE:
%       mat_rs = resize(mat, new_size)
%       mat_rs = resize(mat, new_size, interp_mode)
%
% INPUTS:
%       mat         - matrix to resize
%       new_size    - desired matrix size in elements given by [Nx, Ny] in
%                     2D and [Nx, Ny, Nz] in 3D. Here Nx is the number of
%                     elements in the row direction, Ny is the number of
%                     elements in the column direction, and Nz is the
%                     number of elements in the depth direction.
%
% OPTIONAL INPUTS:
%       interp_mode - interpolation mode used by interp2 and interp3 
%                     (default = '*linear')
%
% OUTPUTS:
%       mat_rs      - resized matrix

% check the inputs for release B.0.2 compatability
if length(varargin{2}) == 1 && nargin >= 3 && length(varargin{3}) == 1

% display warning message
disp('WARNING: input usage deprecated, please see documentation.');
disp('In future releases this usage will no longer be functional.');    

% recursively call resize with the correct inputs
if nargin == 3
    mat_rs = imResample3D(varargin{1}, [varargin{2}, varargin{3}]);
else
    mat_rs = imResample3D(varargin{1}, [varargin{2}, varargin{3}], varargin{4});
end
return

end

% assign the matrix input
mat = varargin{1};

% check for interpolation mode input
if nargin == 2
    interp_mode = '*linear';
elseif nargin ~= 3
    error('incorrect number of inputs');
else
    interp_mode = varargin{3};
end

% check inputs
if ndims(mat) ~= length(varargin{2})
    error('resolution input must have the same number of elements as data dimensions');
end

% perform antialiasing filtering before bilinear or trilinear interpolation
%mat=antialiasing_filter(mat,varargin{2});

switch ndims(mat)
case 2
    % extract the original number of pixels from the size of the matrix
    [Nx_input, Ny_input] = size(mat);

    % extract the desired number of pixels
    Nx_output = varargin{2}(1);
    Ny_output = varargin{2}(2);
        
    % check the size is different to the input size
    if Nx_input ~= Nx_output || Ny_input ~= Ny_output 

        % resize the input matrix to the desired number of pixels
        mat_rs = interp2(0:1/(Ny_input - 1):1, (0:1/(Nx_input - 1):1)', mat, 0:1/(Ny_output - 1):1, (0:1/(Nx_output - 1):1)', interp_mode);

    else
        mat_rs = mat;
    end
case 3

    % extract the original number of pixels from the size of the matrix
    [Nx_input, Ny_input, Nz_input] = size(mat);

    % extract the desired number of pixels
    Nx_output = varargin{2}(1);
    Ny_output = varargin{2}(2); 
    Nz_output = varargin{2}(3);         

    % create normalised plaid grids of current discretisation
    [x_mat, y_mat, z_mat] = ndgrid((0:Nx_input-1)/(Nx_input-1), (0:Ny_input-1)/(Ny_input-1), (0:Nz_input-1)/(Nz_input-1));       

    % create plaid grids of desired discretisation
    [x_mat_interp, y_mat_interp, z_mat_interp] = ndgrid((0:Nx_output-1)/(Nx_output-1), (0:Ny_output-1)/(Ny_output-1), (0:Nz_output-1)/(Nz_output-1));

    % compute interpolation; for a matrix indexed as [M, N, P], the
    % axis variables must be given in the order N, M, P
    mat_rs = interp3(y_mat, x_mat, z_mat, mat, y_mat_interp, x_mat_interp, z_mat_interp, interp_mode); 
case 4
  
    % extract the original number of pixels from the size of the matrix
    [Nx_input, Ny_input, Nz_input, Nc_input] = size(mat);

    % extract the desired number of pixels
    Nx_output = varargin{2}(1);
    Ny_output = varargin{2}(2); 
    Nz_output = varargin{2}(3);  
    Nc_output = Nc_input;

    % create normalised plaid grids of current discretisation
    [x_mat, y_mat, z_mat] = ndgrid((0:Nx_input-1)/(Nx_input-1), (0:Ny_input-1)/(Ny_input-1), (0:Nz_input-1)/(Nz_input-1));       

    % create plaid grids of desired discretisation
    [x_mat_interp, y_mat_interp, z_mat_interp] = ndgrid((0:Nx_output-1)/(Nx_output-1), (0:Ny_output-1)/(Ny_output-1), (0:Nz_output-1)/(Nz_output-1));

    % create output matrix
    mat_rs = zeros(Nx_output,Ny_output,Nz_output,Nc_output);
    
    % loop over all the channel images
    for chn=1:Nc_input
        % compute interpolation; for a matrix indexed as [M, N, P], the
        % axis variables must be given in the order N, M, P
        mat_rs(:,:,:,chn) = interp3(y_mat, x_mat, z_mat, mat(:,:,:,chn), y_mat_interp, x_mat_interp, z_mat_interp, interp_mode); 
    end
    
otherwise
    error('input matrix must be 2,3 or 4 dimensional');
end



function image=antialiasing_filter(image,newsize,method)

if nargin<3
    method='replicate';
end

if (ndims(image)==3)

sc_x =  size(image,1) / newsize(1);
sc_y =  size(image,2) / newsize(2);
sc_z =  size(image,3) / newsize(3);

if ( sc_x <= 1 || sc_y <= 1 || sc_z <= 1)
    return;
end
    
hs_x=fspecial('gaussian',[11,1],sc_x/3.0);
hs_y=fspecial('gaussian',[11,1],sc_y/3.0);
hs_z=fspecial('gaussian',[11,1],sc_z/3.0);

hs_x=reshape(hs_x,[numel(hs_x),1,1]);
image=imfilter(image,hs_x,method);
hs_y=reshape(hs_y,[1,numel(hs_y),1]);
image=imfilter(image,hs_y,method);
hs_z=reshape(hs_z,[1,1,numel(hs_z)]);
image=imfilter(image,hs_z,method);
  
elseif (ndims(image)==2) %#ok<ISMAT>
    
sc_x =  size(image,1) / newsize(1);
sc_y =  size(image,2) / newsize(2);
    
if ( sc_x <= 1 || sc_y <= 1 )
    return;
end
    
hs_x=fspecial('gaussian',[11,1],sc_x/3.0);
hs_y=fspecial('gaussian',[11,1],sc_y/3.0);

hs_x=reshape(hs_x,[numel(hs_x),1]);
image=imfilter(image,hs_x,method);
hs_y=reshape(hs_y,[1,numel(hs_y)]);
image=imfilter(image,hs_y,method);

else
error ('resample::antialiasing dimension error.\n');
end

end



end