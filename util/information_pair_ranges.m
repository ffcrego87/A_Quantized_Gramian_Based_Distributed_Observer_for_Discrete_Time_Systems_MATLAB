function [DO,XO,DM,XM] = information_pair_ranges( ...
    matrices,estimates,p,gammaO,gammaM,rangeFloor,mode)
% Sender-local scales require no maximum-consensus operation.
N = numel(matrices);
n = size(matrices{1},1);
XO = zeros(N,1);
XM = zeros(N,1);
for j = 1:N
    matrix = matrices{j};
    if strcmp(mode,'local'), matrix = 0.5*(matrix+matrix'); end
    message = matrices{j} * estimates((j-1)*n+(1:n));
    XO(j) = gammaO * max(norm(matrix(:),inf),rangeFloor);
    XM(j) = gammaM * max(norm(message(:),inf),rangeFloor);
end
if strcmp(mode,'global')
    XO(:) = max(XO);
    XM(:) = max(XM);
else
    assert(strcmp(mode,'local'),'Unknown range mode.');
end
DO = 2*XO/2^p;
DM = 2*XM/2^p;
end
