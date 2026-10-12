function [Xq,sat_count,code] = uniform_quantizer_sat(X,Delta,Xmax)
%UNIFORM_QUANTIZER_SAT Saturating grid quantizer with an explicit code.
% The grid exponent p is inferred from Delta=2*Xmax/2^p. The signed grid
% index runs from -2^(p-1) to 2^(p-1), inclusive. Adding 2^(p-1) gives
% an unsigned code in 0:2^p, which requires p+1 fixed-length bits.

validateattributes(X,{'numeric'},{'real','finite'});
validateattributes(Delta,{'numeric'},{'real','scalar','positive','finite'});
validateattributes(Xmax,{'numeric'},{'real','scalar','positive','finite'});
sat_mask = abs(X) > Xmax;
sat_count = nnz(sat_mask);

Xsat = min(max(X,-Xmax),Xmax);
p = round(log2((2 * Xmax) / Delta));
signed_index = round(Xsat / Delta);
code = uint64(signed_index + 2^(p-1));
Xq = uniform_quantizer_decode(code,Delta,p);

end
