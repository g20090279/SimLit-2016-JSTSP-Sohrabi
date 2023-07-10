function out = quantizeByCodebook(codeBook, in)
% QUANTIZEBYCODEBOOK quantizes the input to the predefined code word in
% the codebook that has minimum distance in the complex plane.
%
% Input(s):
% - in [scalar, vector, or matrix]: the input data to be quantized 
%
% Output(s):
% - out [in the same dimension of 'in']: the output quantized data
%

codeBook = codeBook(:);

% save the original dimension of input data
[numCol,numRow] = size(in);

in = in(:);
in = repmat(in,1,length(codeBook));

% calculate the distance to each codeword
diff = abs( in - repmat(codeBook.',numCol*numRow,1) );
[~,minInd] = min(diff,[],2);

% the quantized output
out = codeBook(minInd);

% reshape back to original dimension
out = reshape(out, numCol, numRow);

end