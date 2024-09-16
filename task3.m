%% Task 3: Quantization of audio into n number of bits
clearvars
close all

% Read the audio file 'tononoise8mono.wav' and store audio data and sampling rate in y and Fs, respectively
[y, Fs] = audioread('tononoise8mono.wav');

% Store the number of samples to N and generate the corresponding continuous time t
N = length(y) - 1;
t = 0:1/Fs:N/Fs;

% Initialize the number of bits
nbits = 7:-1:1; % for 7 bits down to 1 bit

% Pre-allocate yq matrix to store quantized audio for different bit depths
yq = zeros(length(y), length(nbits) + 1);

% Quantize the audio data into 8-bit signed integer
y8 = round(2^7 * y) / 2^7;  % Quantization for 8 bits
yq(:,1) = y8; % Store y8 as the first column in yq

% Quantize audio data for other bit depths from 7 to 1 bits
for i = 1:length(nbits)
    yq(:,i+1) = round(2^(nbits(i)-1) * y) / 2^(nbits(i)-1); % Quantization for n bits
end

% Plot the quantized audio data for different bit depths
figure(3)
for i = 1:length(nbits) + 1
    subplot(length(nbits) + 1, 1, i);
    plot(t, yq(:,i));
    grid on;
    xlabel('Time (s)');
    ylabel(['yq' num2str(9-i) '[n]']); % Corresponding label for 8 to 1 bits
    if i == 1
        title('Quantized Audio with 8 bits');
    else
        title(['Quantized Audio with ', num2str(nbits(i-1)), ' bits']);
    end
end

% Compute the mean square error (MSE) between the original audio and the quantized versions
mse = zeros(1, length(nbits) + 1);
for i = 1:length(nbits) + 1
    mse(i) = mean((y - yq(:,i)).^2);
end

% Plot the mean square error (MSE)
figure(4)
stem(8:-1:1, mse)
grid on
xlabel('Number of Bits');
ylabel('Mean Square Error (MSE)');
title('Mean Square Error vs Number of Bits');

% Check if MSE is consistent with played back audio (Yes or No)
Response = 'Yes';
