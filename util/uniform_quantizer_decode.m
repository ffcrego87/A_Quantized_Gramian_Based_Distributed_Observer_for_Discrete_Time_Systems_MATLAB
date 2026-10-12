function Xq = uniform_quantizer_decode(code,Delta,p)
%UNIFORM_QUANTIZER_DECODE Decode the p-grid fixed-length representation.
% Valid codewords are the integers 0,...,2^p. A p+1 bit container is used;
% the remaining p+1-bit patterns are unused.

validateattributes(code,{'numeric'},{'integer','nonnegative','<=',2^p});
signed_index = double(code) - 2^(p-1);
Xq = Delta * signed_index;

end
