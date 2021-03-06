function forest = houghForestTrain( data, hs, offsets, varargin )
% Train random forest classifier.
%
% Dimensions:
%  M - number trees
%  F - number features
%  N - number input vectors
%  H - number classes
%  L - number landmarks
%
% USAGE
%  forest = forestTrain( data, hs, [varargin] )
%
% INPUTS
%  data     - [NxF] N length F feature vectors
%  hs       - [Nx1] or {Nx1} target output labels in [1,H]
%  offsets  - [Nx3L]target output landmark offsets
%  varargin - additional params (struct or name/value pairs)
%   .M          - [1] number of trees to train
%   .H          - [max(hs)] number of classes
%   .N1         - [5*N/M] number of data points for training each tree
%   .F1         - [sqrt(F)] number features to sample for each node split
%   .split      - ['gini'] options include 'gini', 'entropy' and 'twoing'
%   .minCount   - [1] minimum number of data points to allow split
%   .minChild   - [1] minimum number of data points allowed at child nodes
%   .maxDepth   - [64] maximum depth of tree
%   .dWts       - [] weights used for sampling and weighing each data point
%   .fWts       - [] weights used for sampling features
%   .discretize - [] optional function mapping structured to class labels
%                    format: [hsClass,hBest] = discretize(hsStructured,H);
%
% OUTPUTS
%  forest   - learned forest model struct array w the following fields
%   .fids     - [Kx1] feature ids for each node
%   .thrs     - [Kx1] threshold corresponding to each fid
%   .child    - [Kx1] index of child for each node
%   .distr    - [KxH] prob distribution at each node
%   .hs       - [Kx1] or {Kx1} most likely label at each node
%   .meanOff  - [Kx3L]or {Kx1} most likely mean landmark offset 
%   .covOff   - [Kx3L^2]or {Kx1} covariance matrix of landmark offsets
%   .count    - [Kx1] number of data points at each node
%   .depth    - [Kx1] depth of each node
%   .gains    - [Kx1] normalized entropy gain for each node 
%               (H_m - w_l*H_l - w_r*H_r) / (prob_reaching_node)

%
% EXAMPLE
%  N=10000; H=5; d=2; [xs0,hs0,xs1,hs1]=demoGenData(N,N,H,d,1,1);
%  xs0=single(xs0); xs1=single(xs1);
%  pTrain={'maxDepth',50,'F1',2,'M',150,'minChild',5};
%  tic, forest=forestTrain(xs0,hs0,pTrain{:}); toc
%  hsPr0 = forestApply(xs0,forest);
%  hsPr1 = forestApply(xs1,forest);
%  e0=mean(hsPr0~=hs0); e1=mean(hsPr1~=hs1);
%  fprintf('errors trn=%f tst=%f\n',e0,e1); figure(1);
%  subplot(2,2,1); visualizeData(xs0,2,hs0);
%  subplot(2,2,2); visualizeData(xs0,2,hsPr0);
%  subplot(2,2,3); visualizeData(xs1,2,hs1);
%  subplot(2,2,4); visualizeData(xs1,2,hsPr1);
%
% See also forestApply, fernsClfTrain
%
% Piotr's Image&Video Toolbox      Version 3.24
% Copyright 2013 Piotr Dollar.  [pdollar-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see external/bsd.txt]

% get additional parameters and fill in remaining parameters
dfs={ 'M',1, 'H',[], 'L',[], 'N1',[], 'F1',[], 'split','gini', 'minCount',1, ...
  'minChild',1, 'maxDepth',64, 'dWts',[], 'fWts',[], 'discretize','', ...
  'regSplit','mse', 'nodeSelectProb',0.5 };
[M,H,L,N1,F1,splitStr,minCount,minChild,maxDepth,dWts,fWts,discretize,regSplitStr,nodeSelectProb] = getPrmDflt(varargin,dfs,1);

[N,F]=size(data); discr=~isempty(discretize); regr=~isempty(offsets); assert(length(hs)==N);
minChild=max(1,minChild); minCount=max([1 minCount minChild]);
if(~regr), nodeSelectProb=0.0; end;
if(regr),  assert(length(offsets)==N);end
if(isempty(L)), L=size(offsets,2)/3; end;
if(isempty(H)), H=max(hs); end; assert(discr || all(hs>0 & hs<=H));
if(isempty(N1)), N1=round(5*N/M); end; N1=min(N,N1);
if(isempty(F1)), F1=round(sqrt(F)); end; F1=min(F,F1);
if(isempty(dWts)), dWts=ones(1,N,'single'); end; dWts=dWts/sum(dWts);
if(isempty(fWts)), fWts=ones(1,F,'single'); end; fWts=fWts/sum(fWts);
split=find(strcmpi(splitStr,{'gini','entropy','twoing'}))-1;
regSplit=find(strcmpi(regSplitStr,{'mse','covariance'}))-1;

if(isempty(split)), error('unknown splitting criteria: %s',splitStr); end
if(isempty(regSplit)), error('unknown regression splitting criteria: %s',regSplit); end
assert ( nodeSelectProb>=0 && nodeSelectProb<=1 );

% make sure data has correct types
if(~isa(data,'single')),           data=single(data);    end
if(~isa(hs,'uint32') && ~discr),     hs=uint32(hs);      end
if(~isa(offsets,'single')),     offsets=single(offsets); end
if(~isa(fWts,'single')),           fWts=single(fWts);    end
if(~isa(dWts,'single')),           dWts=single(dWts);    end

% train M random trees on different subsets of data
prmTree = {H,F1,L,minCount,minChild,maxDepth,fWts,split,discretize,regression,regSplit,nodeSelectProb};
for i=1:M
  if(N==N1), data1=data; hs1=hs; offsets1=offsets; dWts1=dWts; else
    d=wswor(dWts,N1,4); data1=data(d,:); hs1=hs(d); offsets1=offsets(d,:);
    dWts1=dWts(d); dWts1=dWts1/sum(dWts1);
  end
  tree = treeTrain(data1,hs1,offsets1,dWts1,prmTree);
  if(i==1), forest=tree(ones(M,1)); else forest(i)=tree; end
end

end

function tree = treeTrain( data, hs, offsets, dWts, prmTree )

% Train single random tree.
[H,F1,L,minCount,minChild,maxDepth,fWts,split,discretize,regression,regSplit,nodeSelectProb]=deal(prmTree{:});
N=size(data,1); K=2*N-1; discr=~isempty(discretize);
thrs   =zeros(K,1,'single');    distr=zeros(K,H,'single');
fids   =zeros(K,1,'uint32');    gains=zeros(K,1,'single');
meanOff=zeros(K,3*L,'single'); covOff=zeros(K,3*3*L*L,'single');
child=fids; count=fids; depth=fids;
hsn=cell(K,1); dids=cell(K,1); dids{1}=uint32(1:N); k=1; K=2;

while( k < K )
  % get node data and store distribution
  dids1=dids{k}; dids{k}=[]; hs1=hs(dids1); offsets1=offsets(dids1,:);
  n1=length(hs1); count(k)=n1;
  if(discr), [hs1,hsn{k}]=feval(discretize,hs1,H); hs1=uint32(hs1); end
  if(discr), assert(all(hs1>0 & hs1<=H)); end; pure=all(hs1(1)==hs1);
  if(~discr), if(pure), distr(k,hs1(1))=1; hsn{k}=hs1(1); else
      distr(k,:)=histc(hs1,1:H)/n1; [~,hsn{k}]=max(distr(k,:)); end; end
  
  % Non-parametric Parzen estimate of landmark offsets
  if(regression), [meanOff(k,:),covOff(k,:)] = parzenOffset(offsets1); end;
  
  % if pure node or insufficient data don't train split
  if( pure || n1<=minCount || depth(k)>maxDepth ), k=k+1; continue; end
  
  % train split and continue
  fids1=wswor(fWts,F1,4); data1=data(dids1,fids1);
  [~,order1]=sort(data1); order1=uint32(order1-1);
  
  if (rand(1) > nodeSelectProb)
    [fid,thr,gain]=forestFindThr(data1,hs1,dWts(dids1),order1,H,split);
  elseif (regression)
    [fid,thr,gain]=regForestFindThr(data1,offsets1,dWts(dids1),order1,regSplit); 
  else
    error('HoughForestTraining::Node Split Criteria Selection Failed');
  end;
  
  fid=fids1(fid); left=data(dids1,fid)<thr; count0=nnz(left);
  if( gain>1e-10 && count0>=minChild && (n1-count0)>=minChild )
    child(k)=K; fids(k)=fid-1; thrs(k)=thr;
    gains(k)=gain*n1;
    dids{K}=dids1(left); dids{K+1}=dids1(~left);
    depth(K:K+1)=depth(k)+1; K=K+2;
  end; k=k+1;
end

% normalize the computed gain at each node
gains = gains/N;

% create output model struct
K=1:K-1; if(discr), hsn={hsn(K)}; else hsn=[hsn{K}]'; end
tree=struct('fids',fids(K),'thrs',thrs(K),'child',child(K),...
  'distr',distr(K,:),'hs',hsn,'count',count(K),'depth',depth(K),...
  'gains',gains(K), 'meanOff',meanOff(K,:)', 'covOff',covOff(K,:)');
end

function ids = wswor( prob, N, trials )
% Fast weighted sample without replacement. Alternative to:
%  ids=datasample(1:length(prob),N,'weights',prob,'replace',false);
M=length(prob); assert(N<=M); if(N==M), ids=1:N; return; end
if(all(prob(1)==prob)), ids=randperm(M,N); return; end
cumprob=min([0 cumsum(prob)],1); assert(abs(cumprob(end)-1)<.01);
cumprob(end)=1; [~,ids]=histc(rand(N*trials,1),cumprob);
[s,ord]=sort(ids); K(ord)=[1; diff(s)]~=0; ids=ids(K);
if(length(ids)<N), ids=wswor(cumprob,N,trials*2); end
ids=ids(1:N)';
end

function [meanOff,covOff] = parzenOffset(offsets)
meanOff = mean(offsets,1);
covOff  = cov(offsets);
covOff  = covOff(:);
end
