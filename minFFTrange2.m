function [ idx, chunkSamples, chunkSpectrum ] = minFFTrange2( samples, minSize, maxSize, coarseFactor)

    %
    % prepare the first iteration:
    %
    stepFrom    = minSize;
    stepTo      = maxSize;
    stepSize    = (stepTo-stepFrom)/coarseFactor;
    SUMS        = zeros(coarseFactor, 1);
    
    chunksSamples = cell(coarseFactor, 1);
    chunksSpectrum = cell(coarseFactor, 1);
    
    % infinite loop (see break condition...)
    while 1
        
        %
        % Do the process:
        %
        k = 1;
        for chunkSizeFract = stepFrom:stepSize:stepTo
            chunkSize = floor(chunkSizeFract);    
            chunk = samples(1:chunkSize);
            CHUNK = fft(chunk);
            SUMS(k) = sum(abs(CHUNK))/chunkSize / chunkSize;
            chunksSamples{k} = chunk;
            chunksSpectrum{k} = CHUNK;
            k=k+1;
        end
        
        
        % Get the index of the minimum of the fourier transforms:
        [~, mink] = min(SUMS);
        
        % convert the index to index of the samples.
        minIdx = floor(stepFrom+(mink-1)*stepSize);
        
        %
        % break condition:
        %
        if 1 >= stepSize
            %stepSize;
            % plot(samples(1:minIdx*2))
            idx = minIdx;
            chunkSamples = chunksSamples{mink};
            chunkSpectrum = chunksSpectrum{mink};
            break;
        end
        
        %
        % prepare the next iteration:
        %
        
        
        % scan +-1 (current) step in the netx iteration.
        stepFrom = max(minSize, minIdx-stepSize);
        stepTo   = min(maxSize, minIdx+stepSize);
        stepSize = (stepTo-stepFrom)/coarseFactor;
    end

    
end

