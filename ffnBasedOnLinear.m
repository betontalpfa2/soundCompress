% This script reduces the data of a sound.
% This script based on the linearSoundDataCompress.m

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

detuneCoeffWeight = 1.01;
detuneCoeffBias = 0.0001;


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%  CREATE THE NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a new feedforward network with numOfBases hidden layer.
net = feedforwardnet(numOfBases);

% Set the input size equals the size of a chunk.
net.inputs{1}.size = spaceDim;

% Set the output size equals the size of a chunk.
net.layers{2}.size = spaceDim;

net = configure(net,originalSoundMat,originalSoundMat);

net = init(net);

% Set the apriori, linear assumprion for starting.
% (With the detuneCoeffWeight can be started from different places. This
% can be considered as seed.)
net.IW{1} = detuneCoeffWeight.*essentialBases';
net.LW{2, 1} = essentialBases./detuneCoeffWeight;

% % In linear word the biases are 0s:
net.b{1} = ones(numOfBases, 1).*detuneCoeffBias;
net.b{2} = ones(spaceDim, 1).*detuneCoeffBias;

% There are so many weights so trainlm run out of memory. Lets pinc another
% training method.
net.trainFcn = 'traingd';
net.layers{1}.transferFcn = 'elliotsig';
%net.layers{1}.transferFcn = 'purelin';

% Initial training rate parameter for adaptive training. See details below.
net.trainParam.lr  = 1000;
% The adaptive training uses 200 epochs with constant learning parameter,
% then try to adjust it.
net.trainParam.epochs = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%  ADAPTIVE  TRAINING %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adaptive compression.
for k = 1:15
    net.trainParam.lr
    [net, report]  =train(net, originalSoundMat, originalSoundMat);
    if strcmp(report.stop, 'Validation stop.')
        % Validation stops means ca. the exploding gradient problem...
        lastLr = net.trainParam.lr;
        % that the learning rate is too big lets decrease it.
        net.trainParam.lr  = lastLr/3;
    elseif strcmp(report.stop, 'Maximum epoch reached.')
        % Maximum epoch reached means that maybe the learning rate can be
        % inreased.
        lastLr = net.trainParam.lr;
        net.trainParam.lr  = lastLr*3;
    else
        fprintf('ANOTHER STOP!!')
        report.stop
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%  DO THE (DE)COMPRESSION %%%%%%%%%%%%%%%%%%%%%%%

decompressedDataMat = net(originalSoundMat);

decompressedData = reshape(decompressedDataMat, numberOfChunks*spaceDim, 1);
