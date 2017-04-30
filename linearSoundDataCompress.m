% This script reduces the data of a sound.

%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The name of the original sound file.
filename = 'Schubert8.wav';

% This is the dimension of the linear space.
% If there is no compression the number of base vectors will be same.
% The sound will be chunked with this length.
spaceDim = 256;

% The compression is based on that the number of the bases is less than the
% dimension. Because of this not all points in the space could be reached,
% so the compressed sound wont be equal the original. (Lossy compression)
numOfBases = 64;


%%%%%%%%%%%%%%%%%%%%%%%%%%%  PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% open and read the original sound.
originalSound = audioread(filename);

%use only one of the two channel.
originalSound = originalSound(:,1);

% the +-1 computes the ceil...
origLen = size(originalSound, 1);
numberOfChunks = ceil((origLen)/spaceDim);

% fill zeros at the end ot the original Sound in order to have
% numberOfChunks count complete part.
zerosLen = numberOfChunks*spaceDim-origLen;
originalSound = [originalSound; zeros(zerosLen,1)];

% Create a matrix from the chunks.
originalSoundMat = reshape(originalSound, spaceDim, numberOfChunks);

% Compute the covariance of the chunks.
CorrOrigSound = xcorr(originalSound);

firstRow = CorrOrigSound(numberOfChunks*spaceDim:numberOfChunks*spaceDim+spaceDim-1);
CovOrigSound = toeplitz(firstRow);

% compute the eigen vectors and values.
[eigVect, eigVal] = eig(CovOrigSound);

% The eigVect is ordered by the values. The last ones has bigger values.
% keep only the last numOfBases count.
%eigVect(:,1:spaceDim-numOfBases) = zeros(spaceDim,spaceDim-numOfBases);
essentialBases = eigVect(:,spaceDim-numOfBases+1:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%  DO THE COMPRESSION %%%%%%%%%%%%%%%%%%%%%%%%%%%


compressedDataMat = essentialBases'*originalSoundMat;

decompressedDataMat = essentialBases*compressedDataMat;

decompressedData = reshape(decompressedDataMat, numberOfChunks*spaceDim, 1);




