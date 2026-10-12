function [bytes,info] = serialize_information_packet(packet,p,n)
%SERIALIZE_INFORMATION_PACKET Pack an information pair into actual bytes.
% p and n are agreed configuration, not fields transmitted in the packet.
% The first 16 bytes are matrix/vector binary64 ranges in little-endian
% order. Matrix codes (MATLAB column order of triu(true(n))) then vector
% codes follow, each with b=p+1 bits, least-significant bit first. Any last
% byte is padded with zero bits. The decoder rejects unused grid codewords.
% p<=52 keeps every supported codeword exactly representable as a double
% when the existing quantizer reconstructs its signed grid index.

validateattributes(p,{'numeric'},{'scalar','integer','>=',1,'<=',52,'finite'});
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
assert(isstruct(packet) && isscalar(packet),'Expected a scalar packet struct.');
assert(all(isfield(packet,{'matrix_codes','vector_codes','range_matrix','range_vector'})), ...
    'The information packet has missing fields.');
matrixCount = n*(n+1)/2;
validateattributes(packet.matrix_codes,{'numeric'}, ...
    {'real','integer','nonnegative','<=',2^p,'numel',matrixCount});
validateattributes(packet.vector_codes,{'numeric'}, ...
    {'real','integer','nonnegative','<=',2^p,'numel',n});
validateattributes(packet.range_matrix,{'double'},{'real','scalar','positive','finite'});
validateattributes(packet.range_vector,{'double'},{'real','scalar','positive','finite'});

codes = uint64([packet.matrix_codes(:);packet.vector_codes(:)]);
b = p+1;
codeBitsCount = numel(codes)*b;
bitNumber = (0:codeBitsCount-1)';
symbolIndex = floor(bitNumber/b)+1;
bits = uint8(bitget(codes(symbolIndex),mod(bitNumber,b)+1));
paddingBits = mod(-codeBitsCount,8);
bits = [bits;zeros(paddingBits,1,'uint8')];
payload = uint8(sum(double(reshape(bits,8,[])).*(2.^(0:7)'),1))';

rangeBytes = reshape(typecast([packet.range_matrix,packet.range_vector],'uint8'),8,2);
[~,~,endian] = computer;
if endian ~= 'L', rangeBytes = flipud(rangeBytes); end
bytes = [rangeBytes(:);payload];
info = struct('valid_bits',128+codeBitsCount,'padding_bits',paddingBits, ...
    'wire_bytes',numel(bytes),'code_bits_per_scalar',b);
end
