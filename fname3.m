function [b,a,H,h,t,X,y,Y,Yk,MSE,so,H2f,Is,If,Isf]=fname3
          
% De La Salle University
% Electronics and Communications Engineering Department
%
% Course        : LBYCPA4
% SECTION       : EQ1
% Submitted by  : Kenan Banal
% Submitted to  : Dr. Edwin Sybingco
%
% Exercise 3    : Frequency Representations of Signals and Systems
%
% Note: Check the instructions given in canvas

%% Task 1
clc; clearvars;
% Generate the time index for the impulse response
n=0:31;
% Specify the filter coefficients of filter 1 to 4 using the cell array
% b{k} and a{k} where k = 1,2,3, and 4. Also determine the following
% 
% impulse response and store it in h{k}
% Frequency response and store them in H{k}

% Filter 1
b = cell(1, 4); % Cell array to store numerator coefficients
a = cell(1, 4); % Cell array to store denominator coefficients

% H1(z) = 1/8 * (1 + z^(-1)) / (1 - (3/4) * z^(-1))
b{1} = 1/8 * [1 1];        % Numerator coefficients of H1
a{1} = [1 -3/4];           % Denominator coefficients of H1

% H2(z) = 1/8 * (1 - z^(-1)) / (1 + (3/4) * z^(-1))
b{2} = 1/8 * [1 -1];       % Numerator coefficients of H2
a{2} = [1 3/4];            % Denominator coefficients of H2

% H3(z) = 7/32 * (1 - z^(-2)) / (1 + (9/16) * z^(-2))
b{3} = 7/32 * [1 0 -1];    % Numerator coefficients of H3
a{3} = [1 0 9/16];         % Denominator coefficients of H3

% H4(z) = 25/32 * (1 + z^(-2)) / (1 + (9/16) * z^(-2))
b{4} = 25/32 * [1 0 1];    % Numerator coefficients of H4
a{4} = [1 0 9/16];         % Denominator coefficients of H4

% Visualization of figure 1 to 4
H = cell(1, 4);  % Cell array to store frequency response of each filter
w = linspace(0, pi, 512);  % Frequency vector (512 evaluation points)

% Complete the table below by categorizing the filter type as lowpass,
% highpass, bandpass, or bandstop
%figure(1);
%%% 
%%% Plot for Filter 1
%%% 
[H1, w] = freqz(b{1}, a{1}, 512);  % Frequency response of H1
impulse_response1 = filter(b{1}, a{1}, [1 zeros(1, 31)]);  % Impulse response for H1

figure('Name', 'Figure 1: Filter H1', 'Position', [100, 100, 1200, 800]);

% Magnitude Response
subplot(2, 2, 1);
plot(w/pi, abs(H1));
title('Magnitude Response of H1[z]');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude');

% Phase Response
subplot(2, 2, 2);
plot(w/pi, angle(H1));
title('Phase Response of H1[z]');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (radians)');

% Pole-Zero Plot
subplot(2, 2, 3);
zplane(b{1}, a{1});
title('Pole-Zero Plot of H1[z]');

% Impulse Response
subplot(2, 2, 4);
stem(n, impulse_response1, 'filled');
title('Impulse Response of H1[z]');
xlabel('n');
ylabel('Amplitude');
%%% 
%%% Plot for Filter 2
%%% 
[H2, w] = freqz(b{2}, a{2}, 512);  % Frequency response of H2
impulse_response2 = filter(b{2}, a{2}, [1 zeros(1, 31)]);  % Impulse response for H2

figure('Name', 'Figure 2: Filter H2', 'Position', [100, 100, 1200, 800]);

% Magnitude Response
subplot(2, 2, 1);
plot(w/pi, abs(H2));
title('Magnitude Response of H2[z]');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude');

% Phase Response
subplot(2, 2, 2);
plot(w/pi, angle(H2));
title('Phase Response of H2[z]');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (radians)');

% Pole-Zero Plot
subplot(2, 2, 3);
zplane(b{2}, a{2});
title('Pole-Zero Plot of H2[z]');

% Impulse Response
subplot(2, 2, 4);
stem(n, impulse_response2, 'filled');
title('Impulse Response of H2[z]');
xlabel('n');
ylabel('Amplitude');
%%% 
%%% Plot for Filter 3
%%% 
[H3, w] = freqz(b{3}, a{3}, 512);  % Frequency response of H3
impulse_response3 = filter(b{3}, a{3}, [1 zeros(1, 31)]);  % Impulse response for H3

figure('Name', 'Figure 3: Filter H3', 'Position', [100, 100, 1200, 800]);

% Magnitude Response
subplot(2, 2, 1);
plot(w/pi, abs(H3));
title('Magnitude Response of H3[z]');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude');

% Phase Response
subplot(2, 2, 2);
plot(w/pi, angle(H3));
title('Phase Response of H3[z]');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (radians)');

% Pole-Zero Plot
subplot(2, 2, 3);
zplane(b{3}, a{3});
title('Pole-Zero Plot of H3[z]');

% Impulse Response
subplot(2, 2, 4);
stem(n, impulse_response3, 'filled');
title('Impulse Response of H3[z]');
xlabel('n');
ylabel('Amplitude');
%%% 
%%% Plot for Filter 4
%%% 
[H4, w] = freqz(b{4}, a{4}, 512);  % Frequency response of H4
impulse_response4 = filter(b{4}, a{4}, [1 zeros(1, 31)]);  % Impulse response for H4

figure('Name', 'Figure 4: Filter H4', 'Position', [100, 100, 1200, 800]);

% Magnitude Response
subplot(2, 2, 1);
plot(w/pi, abs(H4));
title('Magnitude Response of H4[z]');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude');

% Phase Response
subplot(2, 2, 2);
plot(w/pi, angle(H4));
title('Phase Response of H4[z]');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (radians)');

% Pole-Zero Plot
subplot(2, 2, 3);
zplane(b{4}, a{4});
title('Pole-Zero Plot of H4[z]');

% Impulse Response
subplot(2, 2, 4);
stem(n, impulse_response4, 'filled');
title('Impulse Response of H4[z]');
xlabel('n');
ylabel('Amplitude');

h = [impulse_response1 impulse_response2 impulse_response3 impulse_response3 impulse_response4];

t=table([1;2;3;4],'VariableName',{'FilterNumber'});
t.FilterType(1)={'Highpass'};
t.FilterType(2)={'Lowpass'};
t.FilterType(3)={'Bandstop'};
t.FilterType(4)={'Bandpass'};



%% Task 2: Effect of Number of Samples on Frequency Components
N_values = [25, 50, 75, 100, 200, 1000];  % Different values of N
n_max = 1024;  % Number of points for DTFT (half frequency)
x = cell(1, length(N_values));  % Cell array to store the signals x[n]
X = cell(1, length(N_values));  % Cell array to store the DTFTs X[k]

% Loop over different values of N
for k = 1:length(N_values)
    N = N_values(k);
    n = 0:N-1;  % Time vector for n = 0, 1, ..., N-1
    % Generate the signal x[n]
    x{k} = cos(0.2*pi*n) + cos(0.22*pi*n) + cos(0.6*pi*n);
    
    % Compute the Discrete-Time Fourier Transform (DTFT) with 1024 points
    X{k} = fft(x{k}, n_max);  % Perform FFT and zero-pad to 1024 points
    X{k} = X{k}(1:n_max/2);  % Take only the first half (positive frequencies)
end

% Create figure for magnitude response
figure('Name', 'Figure 5: Magnitude Response for Different N', 'Position', [100, 100, 1200, 800]);

% Plot for each value of N
for k = 1:length(N_values)
    N = N_values(k);
    % Magnitude response for the first half of the DTFT
    freq = linspace(0, 0.5, n_max/2);  % Normalized frequency range (0 to 0.5)
    
    subplot(6, 1, k);  % Create 6 subplots for N = 25, 50, 75, 100, 200, 1000
    plot(freq, abs(X{k}));
    title(['Magnitude Response for N = ', num2str(N)]);
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    ylabel('Magnitude');
    
    % Find peaks with prominence threshold
    [peaks, locs] = findpeaks(abs(X{k}), 'MinPeakProminence', max(abs(X{k}))/2);
    
    % Mark the peaks on the plot
    hold on;
    plot(freq(locs), peaks, 'ro', 'MarkerFaceColor', 'r');
    hold off;
end



%% Task 3 MSE of DTFT and DFT/FFT
% coefficients of the z-transform
%b1=0.0168*[1 0 -2 0 1];
%a1=[1 -2.5333 3.2089 -2.0520 0.6561];
%N=[8 16 32 64];

% Given constants
G = 0.0168;
b = G * [1 0 -2 0 1]; % Numerator coefficients
a = [1 -2.533 3.2089 -2.0520 0.6561]; % Denominator coefficients
N_vals = [8, 16, 32, 64]; % Different values of N
y = cell(1, 4); % Cell array to store impulse responses for each N
Y = cell(1, 4); % Cell array to store DTFT for each N
Yk = cell(1, 4); % Cell array to store DFT for each N
MSE = zeros(1, 4); % Array to store MSE for each N

% Compute the impulse response and frequency responses
for i = 1:length(N_vals)
    N = N_vals(i);
    
    % Compute impulse response using filter (impulse input)
    impulse = [1 zeros(1, N-1)]; % Impulse input of length N
    y{i} = filter(b, a, impulse); % Impulse response
    
    % Compute DTFT using freqz with 1024 points
    [H, w] = freqz(b, a, 1024);
    Y{i} = H; % Store DTFT
    
    % Compute DFT using fft
    Yk{i} = fft(y{i}, 1024); % Zero-pad to 1024 for comparison with DTFT
    
    % Interpolate DFT to match the frequency points of DTFT
    Yk_interp = interp1(linspace(0, pi, N), abs(Yk{i}(1:N)), w(1:N), 'linear', 'extrap');
    
    % Compute MSE between DTFT and DFT
    MSE(i) = mean(abs(abs(Y{i}(1:N)) - Yk_interp).^2); % Mean square error
end

% Visualization (Figure 6) - DTFT overlaying DFT for different N
figure('Name','6','Position', [100, 100, 1200, 800]);
for i = 1:length(N_vals)
    subplot(4, 1, i);
    plot(w/pi, abs(Y{i})); hold on;
    stem(w(1:1024/N_vals(i):1024)/pi, abs(Yk{i}(1:1024/N_vals(i):1024)), 'r'); % Overlay DFT points
    title(['DTFT and DFT overlay for N = ' num2str(N_vals(i))]);
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    ylabel('Magnitude');
    legend('DTFT', 'DFT');
    grid on;
end

% Visualization (Figure 7) - Mean Square Error vs N
figure('Name','7','Position', [100, 100, 1200, 800]);
plot(N_vals, MSE, '-o');
title('Mean Square Error vs N');
xlabel('N');
ylabel('MSE');
grid on;

%% Task 4: time-frequency representation
% Load the audio file
[audio, Fs] = audioread('RollingInTheDeep.wav'); % Fs is the sampling frequency

% Filter coefficients from Task 3
G = 0.0168;
b = G * [1 0 -2 0 1]; % Numerator coefficients
a = [1 -2.533 3.2089 -2.0520 0.6561]; % Denominator coefficients

% Filter the audio signal
so = filter(b, a, audio); % Filtered audio

% Parameters for spectrogram
window = hanning(512); % Hanning window of 512 samples
noverlap = 128; % Overlap of 128 samples
nfft = 512; % 512 frequency points

% Spectrogram of the original audio
figure('Name','8','Position', [100, 100, 1200, 800]);
subplot(2, 1, 1);
spectrogram(audio, window, noverlap, nfft, Fs, 'yaxis');
title('Spectrogram of the Original Audio');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Spectrogram of the filtered audio
subplot(2, 1, 2);
spectrogram(so, window, noverlap, nfft, Fs, 'yaxis');
title('Spectrogram of the Filtered Audio');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

%% Task 5 h2,H2f,I,Is,If,Ifs
% convolution kernel a
h2=0.142*[0 1 2;-1 0 1;-2 -1 0];

% Load the original grayscale image 'lena_gray.tiff'
T = imread('lena_gray.tiff');
T = double(T); % Convert to double for processing

% Define the convolution kernel h2
h2 = 0.142 * [0 1 2; -1 0 1; -2 -1 0];

% Apply convolution
filtered_image = conv2(T, h2, 'same'); % Convolve with the kernel

% Shift amplitude of the filtered image to have a minimum value of 0
filtered_image_shifted = filtered_image - min(filtered_image(:));
filtered_image_uint8 = uint8(filtered_image_shifted); % Convert to uint8

% Store the filtered image in If
If = filtered_image_uint8;

% Compute the 2D Fourier transform of the images and kernel
Is = fftshift(fft2(T)); % Spectrum of the original image
Isf = fftshift(fft2(If)); % Spectrum of the filtered image
H2f = fftshift(fft2(h2, size(T, 1), size(T, 2))); % Spectrum of the convolution kernel

% Convert spectra to logarithmic scale for better visualization
log_Is = log(1 + abs(Is));
log_Isf = log(1 + abs(Isf));
log_H2f = log(1 + abs(H2f));

% Visualization (Figure 9)

figure('Name','9','Position', [100, 100, 1200, 800]);

% Display the original image
subplot(2, 3, 5);
imshow(uint8(T));
title('Original Image');

% Display the filtered image
subplot(2, 3, 6);
imshow(If);
title('Filtered Image');

% Plot the 3D frequency spectrum of the convolution kernel
subplot(2, 3, 1);
mesh(log_H2f);
title('Spectrum of the Filter');
xlabel('Frequency X');
ylabel('Frequency Y');
zlabel('Magnitude (log)');

% Plot the 3D frequency spectrum of the original image
subplot(2, 3, 2);
mesh(log_Is);
title('Spectrum of the Original Image');
xlabel('Frequency X');
ylabel('Frequency Y');
zlabel('Magnitude (log)');

% Plot the 3D frequency spectrum of the filtered image
subplot(2, 3, 3);
mesh(log_Isf);
title('Spectrum of the Filtered Image');
xlabel('Frequency X');
ylabel('Frequency Y');
zlabel('Magnitude (log)');


end



