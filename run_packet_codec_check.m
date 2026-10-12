%% Verify byte serialization, exact reconstruction and malformed packets.
repoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(repoDir,'util'));
rng(173,'twister');
pValues = [1 23 26 28 30 44];
nValues = [2 40];
checked = 0;
for n = nValues
    for p = pValues
        for zeroInput = [false true]
            if zeroInput
                Omega = zeros(n); msg = zeros(n,1);
            else
                Omega = randn(n); Omega = 0.5*(Omega+Omega');
                msg = randn(n,1);
            end
            XO = max(norm(Omega(:),inf),1e-12);
            XM = max(norm(msg,inf),1e-12);
            [canonicalOmega,canonicalMsg,packet] = quantized_information_pair( ...
                Omega,msg,2*XO/2^p,XO,2*XM/2^p,XM);
            [bytes,wire] = serialize_information_packet(packet,p,n);
            [decoded,decodedOmega,decodedMsg,decodedWire] = ...
                deserialize_information_packet(bytes,p,n);
            assert(isequal(decoded,packet),'Packet fields changed in the byte round trip.');
            assert(isequal(decodedOmega,canonicalOmega) && isequal(decodedMsg,canonicalMsg), ...
                'Byte serialization changed a reconstructed numerical value.');
            assert(isequal(wire,decodedWire));
            expectedBits = 128+(n*(n+1)/2+n)*(p+1);
            assert(wire.valid_bits==expectedBits && numel(bytes)==ceil(expectedBits/8));
            assert(wire.padding_bits>=0 && wire.padding_bits<=7);
            checked = checked+1;
        end
    end
end

% Fixed wire fixture: binary64 values 1 and 2 followed by codes 1 and 2
% at b=2. Code bits 1,0,0,1 give the payload byte 9, then zero padding.
fixture = struct('matrix_codes',uint64(1),'vector_codes',uint64(2), ...
    'range_matrix',1.0,'range_vector',2.0);
[fixtureBytes,fixtureWire] = serialize_information_packet(fixture,1,1);
expectedFixture = uint8([0 0 0 0 0 0 240 63 0 0 0 0 0 0 0 64 9])';
assert(isequal(fixtureBytes,expectedFixture),'Wire endianness or code bit order differs.');
assert(fixtureWire.padding_bits==4);

invalidChecks = 0;
for p = pValues
    n = 2;
    packet = struct('matrix_codes',zeros(3,1,'uint64'), ...
        'vector_codes',zeros(2,1,'uint64'),'range_matrix',1.0,'range_vector',1.0);
    [bytes,wire] = serialize_information_packet(packet,p,n);
    badPacket = packet;
    badPacket.matrix_codes(1) = uint64(2^p+1);
    expect_failure(@() serialize_information_packet(badPacket,p,n));
    invalidChecks = invalidChecks+1;

    % Replace the first code (0) by the unused value 2^p+1 in the byte stream.
    badBytes = bytes;
    badBytes(17) = bitor(badBytes(17),uint8(1));
    highByte = 17+floor(p/8);
    badBytes(highByte) = bitset(badBytes(highByte),mod(p,8)+1,1);
    expect_failure(@() deserialize_information_packet(badBytes,p,n));
    invalidChecks = invalidChecks+1;
    expect_failure(@() deserialize_information_packet(bytes(1:end-1),p,n));
    invalidChecks = invalidChecks+1;
    if wire.padding_bits>0
        badBytes = bytes;
        badBytes(end) = bitset(badBytes(end),8,1);
        expect_failure(@() deserialize_information_packet(badBytes,p,n));
        invalidChecks = invalidChecks+1;
    end
end
badFixture = fixture;
badFixture.range_matrix = NaN;
expect_failure(@() serialize_information_packet(badFixture,1,1));
badFixture.range_matrix = 0;
expect_failure(@() serialize_information_packet(badFixture,1,1));
invalidChecks = invalidChecks+2;
fprintf('Packet codec PASS: %i exact reconstruction cases, fixed endian fixture, %i rejected malformed inputs.\n', ...
    checked,invalidChecks);
fprintf('n=40,p=23: valid bits=%i, bytes=%i; n=40,p=25: valid bits=%i, bytes=%i.\n', ...
    128+860*24,(128+860*24)/8,128+860*26,(128+860*26)/8);

function expect_failure(action)
    rejected = false;
    try
        action();
    catch
        rejected = true;
    end
    assert(rejected,'Malformed packet was accepted.');
end
