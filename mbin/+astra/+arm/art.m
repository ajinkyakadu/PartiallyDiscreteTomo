function x = art(W, p, lambda, nIter, x0)

if nargin < 3 || isempty(lambda), lambda = 0.11;              end
if nargin < 4 || isempty(nIter),  nIter = size(W,2)*10;   end
if nargin < 5 || isempty(x0),     x0 = zeros(size(W,2),1); end
    
x = x0;

rowSum = dot(W,W,2);

for i = 1:nIter
    rowIdx = mod(i-1,size(W,1))+1;
    y = (lambda*(p(rowIdx) - W(rowIdx,:)*x)/rowSum(rowIdx)) * W(rowIdx,:);
    x = x + y(:);
%     show(reshape(x, [64,64]));
end
