% This script reduces the data of a sound.
% This script based on linearSoundDataCompress.m
% This script uses multiple sets of bases.
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The name of the original sound file.
filename = 'Schubert8.wav';

% This is the dimension of the linear space.
% If there is no compression the number of base vectors will be same.
% The sound will be chunked to such length.
spaceDim = 512;

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


maxLag = 1023;
VALS = zeros(1, numberOfChunks);
% IDXS = zeros(1, numberOfChunks);
allChunkSamples = cell(1, 5); % the minimum size...
allChunkSpectrum = cell(1, 5); % the minimum size...
minChSize = 512;
maxChSize = 20000;
chunkStartIdx = 1;
i = 1;
IDXS = zeros(1, numberOfChunks);


tmpSound = originalSound;
nextFrom = 1;
fprintf('Start split original sound to chunks...\n');
while minChSize<size(tmpSound, 1) % 1:numberOfChunks
    fprintf('Start cycle: %d\n', i);
    % SUM = ones(1, maxChSize)*1000000; % the lower half will remain zero...
    %chunk = originalSound(chunkStartIdx:chunkStartIdx+minChSize-1);
    %compareTo = originalSound(chunkStartIdx+minChSize/2:chunkStartIdx+maxChSize);
    %chunkCorr = xcorr(compareTo, chunk);
    
    
    % soundsc(originalSound(nextFrom:nextFrom+step), 44100)
    % nextFrom = nextFrom + step;
    maxChSize = min(maxChSize, size(tmpSound, 1));
    [chunkLen, chunkSamples, chunkSpectrum] = minFFTrange2(tmpSound, minChSize, maxChSize, 10);
     
    allChunkSamples{i} = chunkSamples; % chunkSamples;
    allChunkSpectrum{i} = chunkSpectrum;
    tmpSound = tmpSound(chunkLen+1:end);
    
    IDXS(i) = chunkLen;
    % chunkStartIdx = chunkStartIdx+idx;
    i=i+1;
end

if 0<size(tmpSound, 1)
    allChunkSamples{i} = tmpSound;
    allChunkSpectrum{i} = fft(tmpSound);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% RMS EQUALIZER == COMPRESSOR %%%%%%%%%
% len = 2048;
% for i = len:size(originalSound, 1)-len
%     from = i-len/2+1;% max(0, i-)
%     to = from+len-1;
%     RMSval = rms(originalSound(from:to));
%     compSound(i) = originalSound(i)/RMSval;
% end

% testSound = [];


% CHUNKS2 = cell(size(allChunkSamples));

compressedData = cell(size(allChunkSamples));

chunkCorrMat = zeros(size(allChunkSamples, 2));
chunkLengths = zeros(size(allChunkSamples, 2), 1);
fprintf('Start process chunks...\n');
for i = 1:size(allChunkSamples, 2)
    fprintf('Start cycle: %d\n', i);
    
    % Debug / test
    % testSound = [testSound; allChunkSamples{i}];
    
    chunkLengths(i) = size(allChunkSamples{i}, 1);
    
%     for k = 1:size(allChunkSamples, 2)
%         if i<=k
%             chunkCorrMat(i,k) = max(xcorr(allChunkSamples{i}, allChunkSamples{k}));
%             chunkCorrMat(k,i) = chunkCorrMat(i,k);
%         end
%     end
    
    % 
    CHUNK = allChunkSpectrum{i};
    len = max(size(CHUNK));
    CHUNK = CHUNK(1:ceil((len+1)/2));
    
    [values, indexes] = sort(abs(CHUNK), 'descend');
    
    % Lets choose the greatest values.
    % I can see three ways: 
    %  1. Simply choose the 50 biggest one (BUT the length is different so no...)
    %       b) length / 10  (This makes constant bitrate.)
    %  2. Lets pink by RMS. choose the 90% of the chunk rms.
    %  3. Lets choose values above a level.
    
    %
    % 1/b Implemetation: Constant bitrate:
    %
    splitAt = floor(len/8);
    frequencies = indexes(1:splitAt);
    noiseFrequencies = indexes(splitAt+1:end);
    
    
    dataBlock.len = len;
    dataBlock.frequencies = frequencies;
    dataBlock.values = CHUNK(frequencies);
    dataBlock.rms = rms(allChunkSamples{i}); % UNUSED
    dataBlock.noiseLevel = 0; % mean(abs(CHUNK(noiseFrequencies)));
    
    compressedData{i} = dataBlock;
    
end

testSound = [];
fprintf('Start decode blocks...\n');
for i = 1:size(compressedData, 2)
    fprintf('Start cycle: %d\n', i);

    %
    % 1/b Implemetation (see above): Constant bitrate:
    %
    dataBlock = compressedData{i};
    
    CHUNK = ones(dataBlock.len, 1)*dataBlock.noiseLevel;
    CHUNK = CHUNK.*exp(rand(size(CHUNK))*pi*2j);
    %dataBlock.len = len;
    frequencies = dataBlock.frequencies;
    for k=1:max(size(frequencies))
        freq = frequencies(k);
        CHUNK(freq) = dataBlock.values(k);
%         if(freq>1)
%             CHUNK(dataBlock.len-frequencies(k)+2) = conj(dataBlock.values(k));
%         end
    end
    
    halfIdx = ((dataBlock.len+2)/2);
    if rem(halfIdx,1) == 0 % If halfIdx is not a fraction 
        CHUNK(halfIdx) = real(CHUNK(halfIdx))*2;
    end
    CHUNK(ceil(halfIdx):end) = flipud(conj(CHUNK(2:floor(halfIdx))));
    % idx = 19416/2+1
    
    testSound = [testSound; ifft(CHUNK)];
    
end
% sort

ee = originalSound(1:2465000) - testSound(1:2465000);


toc

