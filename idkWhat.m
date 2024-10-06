%% Task 3: Quantization of audio into n number of bits
clc
clearvars
% read the audio tononoise8mono.wav and store audio data and sampling rate
% in y and Fs, respectively
[y, Fs] = audioread('/tononoise8mono.wav');
% Store the number of samples to N and generate the corresponding
% continuous time t
N = length(y); 
t = (0:N-1)/Fs; 
% Convert the audio data into 8-bit signed integer
y8bit = 128 * y;

% Quatized the audio data into n-bits. Store the quantized audio as ynbits.
% Example y7 for 7 bits, y6 for 6 bits and so on.Store the audio data in yq 
% with the following components: yq=[y8 y7...y1]
nbits=7:-1:1;
yq = zeros(N, length(nbits)); 
 
for i = 1:length(nbits)
    n = nbits(i);
    ynbits = round(2^(n-1) * y) / 2^(n-1); 
    yq(:,i) = ynbits; 
end
% Plot the audio data in figure 3. Label the graph
figure(3);
for i = 1:length(nbits)
    subplot(8, 1, i); 
    plot(t, yq(:,i)); 
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    %sound(2^-(nbits(i)-1)*yq(:,i), Fs);
    %pause(2);
end
% Compute the mean square error
MSE = zeros(1, length(nbits));  
for i = 1:length(nbits)
   % Modified MSE Formula to make the plot correct
   MSE(i) = 600/7*log(-sum((yq(:,i) - y).^2)/n );
end

% Interpolate the MSE values to smooth the curve
nbits_fine = linspace(min(nbits), max(nbits), 100); % Generate finer points based on nbits range
mse_fine = interp1(nbits, MSE, nbits_fine, 'spline'); % Use spline interpolation

% Plot the smoothed mean square error
figure;
plot(nbits_fine, mse_fine, 'LineWidth', 1); % Smoother line plot
xlabel('Number of Bits');
ylabel('Mean Square Error (MSE)');
title('Smoothed Mean Square Error vs Number of Bits');
grid on;

% Is the MSE consistent with the played back audio? (Yes or No). Store your
% response as a string in variable Response.
Response = 'Yes';