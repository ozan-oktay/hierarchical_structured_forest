function [IDX,M,ccounts] = meanShiftWeighted(X, radius, weights, rate, maxIter, minCsize, nThreads)
% weighted meanShift clustering algorithm.
%
% Based on the code from Piotr Dollar
%
% USAGE
%  [IDX,M] = meanShift(X, radius, [rate], [maxIter], [minCsize], [blur] )
%
% INPUTS
%  X           - column vector of data - N vectors of dim p (X is Nxp)
%  radius      - the bandwidth (radius of the window)
%  weights     - column vector of data weights - (Nx1)
%  rate        - [] gradient descent proportionality factor in (0,1]
%  maxIter     - [] maximum number of iterations
%  minCsize    - [] min cluster size (smaller clusters get eliminated)
%
% OUTPUTS
%  IDX         - cluster membership [see kmeans2.m]
%  M           - cluster means [see above]
%
% EXAMPLE
%
% See also MEANSHIFTIM, DEMOCLUSTER

if( nargin<4 ); rate =.2; end
if( nargin<5 ); maxIter =100; end
if( nargin<6 ); minCsize = 1; end
if( nargin<7 ); nThreads = 1; end
if( rate<=0 || rate>1 ); error('rate must be between 0 and 1'); end

% c code does the work  (meanShift1 requires X')
[IDX,meansFinal] = meanShiftWeighted1(X',radius,rate,maxIter,weights,nThreads);
meansFinal = meansFinal';

% calculate final cluster means per cluster
p = size(X,2);  k = max(IDX); M = zeros(k,p); 
for i=1:k 
  idx_m  = (IDX==i);
  M(i,:) = sum(meansFinal(idx_m,:).*repmat(weights(idx_m),1,p)) / sum(weights(idx_m));
end

% sort clusters [largest first] and remove all smaller then minCsize
ccounts = zeros(1,k); 
for i=1:k
  ccounts(i) = sum( weights(IDX==i) ); 
end

[ccounts,order] = sort( -ccounts ); ccounts = -ccounts; M = M(order,:);
IDX2 = IDX;  for i=1:k; IDX2(IDX==order(i))=i; end; IDX = IDX2;
[v,loc] = min( ccounts>=minCsize );
if( v==0 ); M( loc:end, : ) = []; IDX( IDX>=loc ) = -1; end