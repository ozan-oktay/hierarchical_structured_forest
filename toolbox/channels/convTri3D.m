function J = convTri3D( I, r, convdim, s)

if( nargin<3 ), convdim=ones(1,min(3,ndims(I))); end
if( nargin<4 ), s=1; end
if( isempty(I) || (r==0 && s==1) ), J = I; return; end

if (ndims(I)==3 || ndims(I)==4)
 
  if(r<=1), p=12/r/(r+2)-2; f=[1 p 1]/(2+p); r=1;
  else f=[1:r r+1 r:-1:1]/(r+1)^2; end
  J = padarray(I,(convdim*r),'symmetric','both');

  hs_x=reshape(f,[numel(f),1,1]);
  hs_y=reshape(f,[1,numel(f),1]);
  hs_z=reshape(f,[1,1,numel(f)]);

  if(convdim(1)), J = convn(J,hs_x,'valid'); end;
  if(convdim(2)), J = convn(J,hs_y,'valid'); end;
  if(convdim(3)), J = convn(J,hs_z,'valid'); end;
  
  if(s>1), t=floor(s/2)+1;
    if(convdim(1)), J=J(t:s:end-s+t,:,:,:); end;
    if(convdim(2)), J=J(:,t:s:end-s+t,:,:); end;
    if(convdim(3)), J=J(:,:,t:s:end-s+t,:); end;
  end

elseif (ndims(I)==2) %#ok<ISMAT>
          
  if(r<=1), p=12/r/(r+2)-2; f=[1 p 1]/(2+p); r=1;
  else f=[1:r r+1 r:-1:1]/(r+1)^2; end
  J = padarray(I,[r r],'symmetric','both');     
     
  hs_x=reshape(f,[numel(f),1,1]);
  hs_y=reshape(f,[1,numel(f),1]);

  J = convn(J,hs_x,'valid');
  J = convn(J,hs_y,'valid');
  
  if(s>1), t=floor(s/2)+1; J=J(t:s:end-s+t,t:s:end-s+t,:); end

else
  error ('convTri::dimension error');
          
end

