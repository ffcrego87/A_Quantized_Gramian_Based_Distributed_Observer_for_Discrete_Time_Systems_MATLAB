function [packet,Omegaq,msgq,info] = deserialize_information_packet(bytes,p,n)
%DESERIALIZE_INFORMATION_PACKET Read the documented little-endian wire format.
% p and n are agreed configuration. Exactly one packet must be supplied.
validateattributes(p,{'numeric'},{'scalar','integer','>=',1,'<=',52,'finite'});
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
matrixCount = n*(n+1)/2;
b = p+1;
codeBitsCount = (matrixCount+n)*b;
expectedBytes = 16+ceil(codeBitsCount/8);
validateattributes(bytes,{'uint8'},{'vector','numel',expectedBytes});
bytes = bytes(:);

rangeBytes = reshape(bytes(1:16),8,2);
[~,~,endian] = computer;
if endian ~= 'L', rangeBytes = flipud(rangeBytes); end
ranges = typecast(rangeBytes(:),'double');
validateattributes(ranges,{'double'},{'real','positive','finite','numel',2});

payload = bytes(17:end);
bitNumber = (0:8*numel(payload)-1)';
bits = bitget(payload(floor(bitNumber/8)+1),mod(bitNumber,8)+1);
assert(all(bits(codeBitsCount+1:end)==0),'Nonzero packet padding bits.');
codeBits = reshape(bits(1:codeBitsCount),b,[]);
codes = uint64(sum(double(codeBits).*(2.^(0:b-1)'),1))';
validateattributes(codes,{'uint64'},{'<=',2^p});
packet = struct('matrix_codes',codes(1:matrixCount), ...
    'vector_codes',codes(matrixCount+1:end), ...
    'range_matrix',ranges(1),'range_vector',ranges(2));

mask = triu(true(n));
Omegaq = zeros(n);
Omegaq(mask) = uniform_quantizer_decode(packet.matrix_codes,2*ranges(1)/2^p,p);
Omegaq = Omegaq+triu(Omegaq,1)';
msgq = uniform_quantizer_decode(packet.vector_codes,2*ranges(2)/2^p,p);
info = struct('valid_bits',128+codeBitsCount, ...
    'padding_bits',8*numel(payload)-codeBitsCount, ...
    'wire_bytes',expectedBytes,'code_bits_per_scalar',b);
end
