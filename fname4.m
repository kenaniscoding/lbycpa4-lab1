function [w,W,b1,H1,b2,H2,H3,y,S]=fname4
% De La Salle University
% Electronics and Computer Engineering Department
%
% Course        : LBYEC4A/LBYCPA4
% SECTION       : EQ1
% Submitted by  : Kenan Banal
% Submitted to  : Dr. Edwin Sybingco
%
% Exercise 4    : FIR Filter
% Note: Check the instructions given in canvas
%% Task 1
% Specify the order of the window
M = 90;
% Generate the time index
n = 0:M;

% Generate the windows and store them in the matrix w
w = zeros(M+1, 5); % Initialize a matrix with 91 rows and 5 columns

% 1st column: Boxcar window
w(:, 1) = boxcar(M+1);

% 2nd column: Hanning window
w(:, 2) = hanning(M+1);

% 3rd column: Hamming window
w(:, 3) = hamming(M+1);

% 4th column: Blackman window
w(:, 4) = blackman(M+1);

% 5th column: Kaiser window with beta = 3.3953
beta = 3.3953;
w(:, 5) = kaiser(M+1, beta);

% Generate the frequency response for each window and store them in matrix W
numFreqPoints = 512;
W = zeros(numFreqPoints, 5); % Initialize frequency response matrix

for i = 1:5
    % Calculate the frequency response of each window
    [H, f] = freqz(w(:, i), 1, numFreqPoints, 'half');
    W(:, i) = H; % Store the frequency response magnitude
end

% Plot Figure 1: Continuous plot of all the windows
figure;
plot(n, w);
title('Figure 1: Window Functions');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Boxcar', 'Hanning', 'Hamming', 'Blackman', 'Kaiser');
grid on

% Plot Figure 2: Magnitude and Phase Spectrum
figure;

% Magnitude Spectrum in dB
subplot(2, 1, 1);
plot(f, 20*log10(abs(W)));
title('Magnitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Boxcar', 'Hanning', 'Hamming', 'Blackman', 'Kaiser');
grid on

% Phase Spectrum
subplot(2, 1, 2);
plot(f, angle(W));
title('Phase Spectrum');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
legend('Boxcar', 'Hanning', 'Hamming', 'Blackman', 'Kaiser');
grid on


%% Task 2
% Specify parameters for the fir2 function
M = 90;  % Order of the filter
n = 0:M; % Time index

% Define the frequency vector and desired magnitude response for the ideal filter
% This frequency response can be adjusted as per specific requirements
freq_points = [0, 0.2, 0.2, 0.5, 0.5, 0.8, 0.8, 1.0];
desired_mag = [0.2, 0.2, 1.0, 1.0, 0.6, 0.6, 0.8, 0.8];

% Initialize matrix b1 to store the filter coefficients of each windowed filter
b1 = zeros(5, M+1);

% Design five filters using the fir2 function and each window from Task 1
for i = 1:5
    b1(i, :) = fir2(M, freq_points, desired_mag) .* w(:, i).';
end

% Determine the frequency response of each filter and store in H1
numFreqPoints = 512;
H1 = zeros(numFreqPoints, 5); % Initialize matrix to store frequency response

for i = 1:5
    % Calculate the frequency response of each filter
    [H, f] = freqz(b1(i, :), 1, numFreqPoints, 'half');
    H1(:, i) = H; % Store the frequency response
end
name= ["boxcar" "hanning" "hamming" "blackman" "kaiser"];
% Plot Figure 3: Impulse response of the five filters
figure('Position', [100, 100, 1200, 800]);
for i = 1:5
    subplot(5, 1, i);
    stem(n, b1(i, :));
    title(['Figure 3: ', name(i)]);
    xlabel('Sample Index');
    ylabel('Amplitude');
    grid on
end

% Plot Figure 4: Magnitude and Phase Spectrum
figure;
f6 = linspace(0, 1, 512); % Frequency vector for plotting
% Magnitude Spectrum in dB
subplot(2, 1, 1);
plot(f6, log(abs(H1)));
title('Magnitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude ');
legend('boxcar', 'hanning', 'hamming', 'blackman', 'kaiser');
grid on
% Phase Spectrum
subplot(2, 1, 2);
plot(f6, angle(H1));
title('Phase Spectrum');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
legend('boxcar', 'hanning', 'hamming', 'blackman', 'kaiser');
grid on

%% Task 3

% Define filter specifications
delta1 = 0.0575; % Passband deviation
delta2_lowpass = 0.0009; % Stopband deviation for lowpass
delta2_bandpass = 0.0009; % Stopband deviation for bandpass
delta2_highpass = 0.0016; % Stopband deviation for highpass

% Convert gain deviations to decibel values
Rp = 20 * log10((1 + delta1) / (1 - delta1));
Rs_lowpass = -20 * log10(delta2_lowpass);
Rs_bandpass = -20 * log10(delta2_bandpass);
Rs_highpass = -20 * log10(delta2_highpass);

% Define passband and stopband edges in normalized frequency (0 to 1)
wp_lowpass = 0.25 * pi;
ws_lowpass = 0.3 * pi;

wp_bandpass = [0.3 * pi, 0.75 * pi];
ws_bandpass = [0.25 * pi, 0.8 * pi];

wp_highpass = 0.8 * pi;
ws_highpass = 0.75 * pi;

% Normalize frequencies for firpmord
wp_lowpass = wp_lowpass / pi;
ws_lowpass = ws_lowpass / pi;

wp_bandpass = wp_bandpass / pi;
ws_bandpass = ws_bandpass / pi;

wp_highpass = wp_highpass / pi;
ws_highpass = ws_highpass / pi;

% Calculate filter order and design filters using firpmord and firpm
% Initialize variables
N = zeros(1, 3); % Filter orders
b2 = {}; % Cell array to store filter coefficients for different types

% Lowpass filter design
[N(1), fo, ao, w] = firpmord([wp_lowpass, ws_lowpass], [1, 0], [delta1, delta2_lowpass]);
b2{1} = firpm(N(1), fo, ao, w);

% Bandpass filter design
[N(2), fo, ao, w] = firpmord([ws_bandpass(1), wp_bandpass(1), wp_bandpass(2), ws_bandpass(2)], ...
    [0, 1, 0], [delta2_bandpass, delta1, delta2_bandpass]);
b2{2} = firpm(N(2), fo, ao, w);

% Highpass filter design
[N(3), fo, ao, w] = firpmord([ws_highpass, wp_highpass], [0, 1], [delta2_highpass, delta1]);
b2{3} = firpm(N(3), fo, ao, w);

% Calculate the frequency response for each filter
numFreqPoints = 512;
H2 = zeros(numFreqPoints, 3); % Initialize matrix to store frequency response
f = linspace(0, 1, numFreqPoints); % Frequency vector for plotting

for i = 1:3
    [H, ~] = freqz(b2{i}, 1, numFreqPoints, 'half');
    H2(:, i) = H; % Store the frequency response
end

%% Task 3

% Specifications for the filters
delta1 = 0.0575;   % Passband ripple
delta2_low = 0.0009;  % Lowpass stopband deviation
delta2_band = 0.0009; % Bandpass stopband deviation
delta2_high = 0.0016; % Highpass stopband deviation

% Convert passband/stopband deviations to decibels
Rp = 20 * log10((1 + delta1) / (1 - delta1));  % Passband ripple in dB
Rs_low = -20 * log10(delta2_low);              % Lowpass stopband attenuation in dB
Rs_band = -20 * log10(delta2_band);            % Bandpass stopband attenuation in dB
Rs_high = -20 * log10(delta2_high);            % Highpass stopband attenuation in dB

% Define passband and stopband edge frequencies for each filter
wp_low = 0.25 * pi; ws_low = 0.3 * pi;         % Lowpass filter frequencies
wp_band = [0.3 * pi, 0.75 * pi]; ws_band = [0.25 * pi, 0.8 * pi]; % Bandpass
wp_high = 0.8 * pi; ws_high = 0.75 * pi;       % Highpass filter frequencies

% Determine the required order and design filters using firpm
N = zeros(1, 3);       % Initialize vector to store filter orders
b2 = {};               % Initialize cell array to store filter coefficients

% Lowpass filter
[N(1), fo, ao, w] = firpmord([wp_low, ws_low] / pi, [1, 0], [delta1, delta2_low]);
b2{1} = firpm(N(1), fo, ao, w);

% Bandpass filter
[N(2), fo, ao, w] = firpmord([ws_band(1), wp_band(1), wp_band(2), ws_band(2)] / pi, [0, 1, 0], [delta2_band, delta1, delta2_band]);
b2{2} = firpm(N(2), fo, ao, w);

% Highpass filter
[N(3), fo, ao, w] = firpmord([ws_high, wp_high] / pi, [0, 1], [delta2_high, delta1]);
b2{3} = firpm(N(3), fo, ao, w);

% Calculate frequency responses of the filters
numFreqPoints = 512;
H2 = zeros(numFreqPoints, 3);  % Matrix to store frequency responses

for i = 1:3
    [H, f] = freqz(b2{i}, 1, numFreqPoints, 'half');
    H2(:, i) = H; % Store frequency response
end

titlesF = ["Lowpass", "Bandpass", "Highpass"];
% Figure 5: Impulse Response
figure;
for i = 1:3
    subplot(3, 1, i);
    stem(0:N(i), b2{i});
    title('Figure 5: ', titlesF(i));
    xlabel('Sample Index');
    ylabel('Amplitude');
end

% Figure 6: Magnitude and Phase Spectrum
figure;
f = linspace(0, 1, numFreqPoints); % Frequency vector for plotting
% Magnitude Spectrum in dB
subplot(2, 1, 1);
plot(f, (abs(H2)));
title('Magnitude Spectrum ');
xlabel('Frequency (Hz)');
ylabel('Magnitude ');
legend('Lowpass', 'Bandpass', 'Highpass');
grid on

% Phase Spectrum
subplot(2, 1, 2);
plot(f, angle(H2));
title('Phase Spectrum');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
legend('Lowpass', 'Bandpass', 'Highpass');
grid on


%% Task 4
alpha =[1 1 1;...
        1 0.5 1;...
        0.5 1 0.5;...
        1 0.6 0.2];
% Initialize H3 to store the frequency responses for each setting
numFreqPoints = 512;
H3 = zeros(numFreqPoints, 4);

% Calculate the combined frequency response for each setting
for k = 1:4
    % Compute the weighted sum of the individual filter responses
    H3(:, k) = alpha(k, 1) * H2(:, 1) + alpha(k, 2) * H2(:, 2) + alpha(k, 3) * H2(:, 3);
end

%% Visualization

% Plot Figure 7: Magnitude and Phase Response for each setting
figure;

f7 = linspace(0, 1, 512); % Frequency vector for plotting

for k = 1:4
    % Magnitude Response in dB for each setting
    subplot(2, 4, k);
    plot(f7, abs(H3(:, k)));
    title(['Figure 7: Magnitude Response (Setting ', num2str(k), ')']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on
    % Phase Response for each setting
    subplot(2, 4, k+4);
    plot(f7, angle(H3(:, k)));
    title(['Figure 7: Phase Response (Setting ', num2str(k), ')']);
    xlabel('Frequency (Hz)');
    ylabel('Phase (radians)');
    grid on
end


%% Task 5

% Read the audio file and store it in the first column of y
[y(:,1), Fs] = audioread('RollingInTheDeep.wav');  % Original audio

% Apply the three filters from Task 4 (Setting 1: alpha = [1 1 1])
y(:,2) = filter(b2{1}, 1, y(:,1));  % Lowpass filter output
y(:,3) = filter(b2{2}, 1, y(:,1));  % Bandpass filter output
y(:,4) = filter(b2{3}, 1, y(:,1));  % Highpass filter output

% Spectrogram settings
window = kaiser(256, 3.3953);  % Kaiser window with beta = 3.3953
noverlap = 128;                % Overlap of 128 samples
nfft = 512;                    % 512-point FFT for half-frequency response

% Initialize the matrix S to store spectrograms of each channel in y
[S, ~, ~] = spectrogram(y(:,1), window, noverlap, nfft, Fs);  % Pre-allocate size based on 1st spectrogram
S = zeros(size(S, 1), size(S, 2), 4);  % S will be 257 x 1032 x 4

% Calculate the spectrogram for each channel
for k = 1:4
    [S(:,:,k), f, t] = spectrogram(y(:,k), window, noverlap, nfft, Fs);
end

% Plot Figure 8: Spectrograms in a 3D scatter plot
figure('Position', [100, 100, 1200, 800]);
titles = ["Original", "Lowpass", "Bandpass", "Highpass"];

for k = 1:4
    % Prepare data for 3D scatter plot
    [T, F] = meshgrid(t, f);                     % Time-frequency grid
    Magnitude = 20*log10(abs(S(:,:,k)));         % Magnitude in dB
    
    % Map the Magnitude to a color based on the current colormap
    colormap jet;                                % Use the 'jet' colormap
    clim([-120, 0]);                            % Set color axis limits (adjust as needed)
    colors = colormap;                           % Get current colormap
    colorIndices = round(rescale(Magnitude(:), 1, size(colors, 1))); % Scale magnitudes to colormap index
    
    % Scatter plot for each spectrogram
    subplot(2, 2, k);
    scatter3(T(:), F(:), Magnitude(:), 1, colors(colorIndices, :), 'filled');  % Plot with color-mapped dots
    title([titles(k), ' Spectrogram']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    zlabel('Magnitude');
    grid on
    colorbar;
end




end
