clear all
clc

% Read the input audio file
[inputAudio, Fs] = audioread('mono.wav');

% Define the window size and overlap
windowSize = 2; % Can be modified to evaluate different sizes
overlap = windowSize/2; % 50% overlap

% Pre-process the audio signal
numSamples = length(inputAudio);
numWindows = floor((numSamples - overlap) / (windowSize - overlap));

% Initialize encoded matrix
encodedAudio = zeros(windowSize, numWindows);

% Quantization level (adjust this to change the bitrate)
quantizationLevels = 64; % Example value, adjust as needed

% Apply windowing and MDCT
for i = 1:numWindows
    startIdx = (i - 1) * (windowSize - overlap) + 1;
    endIdx = startIdx + windowSize - 1;
    
    % Extract the windowed segment
    windowedSegment = inputAudio(startIdx:endIdx);
    
    % Apply MDCT
    encodedSegment = mdct(windowedSegment);
    
    % Quantize the MDCT coefficients
    maxVal = max(abs(encodedSegment));
    quantizedSegment = round(encodedSegment / maxVal * (quantizationLevels / 2)) * (maxVal / (quantizationLevels / 2));
    
    % Store the quantized segment
    encodedAudio(:, i) = quantizedSegment;
end

% To decode the audio, apply inverse MDCT (IMDCT) and overlap-add technique
decodedAudio = zeros(numSamples, 1);

for i = 1:numWindows
    startIdx = (i - 1) * (windowSize - overlap) + 1;
    endIdx = startIdx + windowSize - 1;
    
    % Extract the encoded segment
    quantizedSegment = encodedAudio(:, i);
    
    % Apply IMDCT
    decodedSegment = imdct(quantizedSegment);
    
    % Overlap and add to the decoded audio
    decodedAudio(startIdx:endIdx) = decodedAudio(startIdx:endIdx) + decodedSegment;
end

% Normalize the decoded audio to prevent clipping
decodedAudio = decodedAudio / max(abs(decodedAudio));

% Save the decoded audio to a file
audiowrite('output.wav', decodedAudio, Fs);

% Calculate SNR
signalPower = sum(inputAudio .^ 2) / numSamples;
noisePower = sum((inputAudio - decodedAudio) .^ 2) / numSamples;

% Ensure noisePower is not zero to avoid NaN SNR
if noisePower == 0
    SNR = Inf;
else
    SNR = 10 * log10(signalPower / noisePower);
end

% Calculate compression ratio
audioFile = dir('mono.wav');
originalSize = audioFile.bytes;
%originalSize = numSamples * 16 / 1024; % Original size in kb (assuming 16-bit audio)
encodedSize = numWindows * windowSize * log2(quantizationLevels) /8 /1024; % Encoded size in kb
compressionRatio = originalSize / (encodedSize*1024);

% Calculate size of encoded audio in kb
encodedAudioSize = encodedSize;

% Display the results
fprintf('SNR: %.2f dB\n', SNR);
fprintf('Compression Ratio: %.2f\n', compressionRatio);
fprintf('Size of Encoded Audio: %.2f kb\n', encodedAudioSize);

% Plot the original and reconstructed signals
figure;
subplot(2,1,1);
plot(inputAudio);
title('Original Audio Signal');
subplot(2,1,2);
plot(decodedAudio);
title('Reconstructed Audio Signal');

function X = mdct(x)
    N = length(x);
    w = sin(pi/(2*N) * ((1:N) - 0.5))';
    x = x .* w;
    X = zeros(N, 1);
    for k = 0:N-1
        X(k+1) = sum(x .* cos(pi/N * (k + 0.5) * ((1:N) - 0.5))');
    end
end

function x = imdct(X)
    N = length(X);
    w = sin(pi/(2*N) * ((1:N) - 0.5))';
    x = zeros(N, 1);
    for n = 0:N-1
        x(n+1) = sum(X .* cos(pi/N * ((1:N) - 0.5) * (n + 0.5))');
    end
    x = x .* w;
end
