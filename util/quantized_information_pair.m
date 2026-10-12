function [Omegaq,msgq,packet,diagnostic] = quantized_information_pair( ...
    Omega,msg,DeltaOmega,XmaxOmega,DeltaMsg,XmaxMsg)
% Encode one symmetric triangle and the vector; transmit two binary64 ranges.
% The decoder uses only packet fields and the agreed grid exponent p.
Omega = 0.5 * (Omega + Omega');
n = size(Omega,1);
mask = triu(true(n));
[canonical,satFull] = uniform_quantizer_sat(Omega,DeltaOmega,XmaxOmega);
[~,satTriangle,matrixCodes] = uniform_quantizer_sat(Omega(mask),DeltaOmega,XmaxOmega);
[canonicalMsg,satMsg,vectorCodes] = uniform_quantizer_sat(msg,DeltaMsg,XmaxMsg);
p = round(log2(2 * XmaxOmega / DeltaOmega));
assert(p == round(log2(2 * XmaxMsg / DeltaMsg)));
packet = struct('matrix_codes',matrixCodes,'vector_codes',vectorCodes, ...
    'range_matrix',double(XmaxOmega),'range_vector',double(XmaxMsg));
[wireBytes,wireInfo] = serialize_information_packet(packet,p,n);
[decoded,Omegaq,msgq] = deserialize_information_packet(wireBytes,p,n);
assert(isequal(packet,decoded),'Byte packet reconstruction failed.');
verified = isequal(Omegaq,canonical) && isequal(msgq,canonicalMsg);
assert(verified,'Packed reconstruction differs from canonical full quantization.');
diagnostic = struct('full_matrix_saturations',satFull, ...
    'matrix_code_saturations',satTriangle,'vector_saturations',satMsg, ...
    'roundtrip_verified',verified,'payload_bits',wireInfo.valid_bits, ...
    'wire_bytes',wireInfo.wire_bytes,'padding_bits',wireInfo.padding_bits);
end
