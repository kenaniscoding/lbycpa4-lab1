function [b1,a1,H1,Grp1,h1,b2,a2,H2,Grp2,h2,y,S]=fname5
% De La Salle University
% Department of Electronics and Computer Engineering 
%
% Course        : LBYCPA4
% SECTION       : EQ1
% Submitted by  : Kenan Banal
% Submitted to  : Dr.Edwin Sybingco
%
% Exercise 5    : IIR Filter Design

% Load the audio
[s,Fs]=audioread('RollingInTheDeep.wav');
Rp=0.5;
Rs=30;


% Define the sampling frequency
Fs = 44100;

% Define passband and stopband frequencies for each filter
F_pass = [4100, 4250, 8950, 13400, 17000];
F_stop = [4500, 3900, 8350, 13000, 17400];
F_pass2 = [NaN, 8750, 13250, 16800, NaN];
F_stop2 = [NaN, 9350, 13650, 17500, NaN];


% Initialize cell arrays to store coefficients and responses
b1 = cell(5, 1); % Numerator coefficients
a1 = cell(5, 1); % Denominator coefficients
H1 = cell(5, 1); % Frequency response
Grp1 = cell(5, 1); % Group delay
h1 = cell(5, 1); % Impulse response

% Loop through each filter and design it
for i = 1:5
    if isnan(F_pass2(i))
        % Lowpass or Highpass Filter
        if i == 1
            [N, Wn] = cheb1ord(F_pass(i)/(Fs/2), F_stop(i)/(Fs/2), Rp, Rs);
            [b1{i}, a1{i}] = cheby1(N, Rp, Wn, 'low');
        elseif i == 5
            [N, Wn] = cheb1ord(F_pass(i)/(Fs/2), F_stop(i)/(Fs/2), Rp, Rs);
            [b1{i}, a1{i}] = cheby1(N, Rp, Wn, 'high');
        end
    else
        % Bandpass Filter
        [N, Wn] = cheb1ord([F_pass(i) F_pass2(i)]/(Fs/2), ...
                           [F_stop(i) F_stop2(i)]/(Fs/2), Rp, Rs);
        [b1{i}, a1{i}] = cheby1(N, Rp, Wn, 'bandpass');
    end
    
    % Frequency response
    [H1{i}, w] = freqz(b1{i}, a1{i}, 1024, Fs);
    
    % Group delay
    Grp1{i} = grpdelay(b1{i}, a1{i}, 1024, Fs);
    
    % Impulse response
    h1{i} = impz(b1{i}, a1{i}, 120);
end
% Define colors for each filter
colors = lines(5); % Generates a set of distinct colors for plotting

% Combined Magnitude Spectrum Plot (Linear Scale)
figure(1);
hold on;
for i = 1:5
    plot(w, abs(H1{i}), 'Color', colors(i,:), 'DisplayName', ['Filter ', num2str(i)]);
end
hold off;
title('Combined Magnitude Response of Filters (Linear Scale)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (Linear Scale)');
ylim([0, 1]); % Set y-axis range from 0 to 1
legend show;
grid on;

% Combined Phase Spectrum Plot with y-axis range -4 to 4
figure(2);
hold on;
for i = 1:5
    plot(w, angle(H1{i}), 'Color', colors(i,:), 'DisplayName', ['Filter ', num2str(i)]);
end
hold off;
title('Combined Phase Response of Filters');
xlabel('Frequency (Hz)');
ylabel('Phase (Radians)');
ylim([-4, 4]); % Set y-axis range from -4 to 4
legend show;
grid on;

% Combined Group Delay Plot for Chebyshev Type 1 Filters
figure(3);
hold on;
for i = 1:5
    plot(w, Grp1{i}, 'Color', colors(i,:), 'DisplayName', ['Filter ', num2str(i)]);
end
hold off;
title('Combined Group Delay of Chebyshev Type 1 Filters');
xlabel('Frequency (Hz)');
ylabel('Group Delay (samples)');
legend show;
grid on;


% Initialize cell arrays to store coefficients and responses
b2 = cell(5, 1); % Numerator coefficients
a2 = cell(5, 1); % Denominator coefficients
H2 = cell(5, 1); % Frequency response
Grp2 = cell(5, 1); % Group delay
h2 = cell(5, 1); % Impulse response

% Loop through each filter and design it
for i = 1:5
    if isnan(F_pass2(i))
        % Lowpass or Highpass Filter
        if i == 1
            [N, Wn] = ellipord(F_pass(i)/(Fs/2), F_stop(i)/(Fs/2), Rp, Rs);
            [b2{i}, a2{i}] = ellip(N, Rp, Rs, Wn, 'low');
        elseif i == 5
            [N, Wn] = ellipord(F_pass(i)/(Fs/2), F_stop(i)/(Fs/2), Rp, Rs);
            [b2{i}, a2{i}] = ellip(N, Rp, Rs, Wn, 'high');
        end
    else
        % Bandpass Filter
        [N, Wn] = ellipord([F_pass(i) F_pass2(i)]/(Fs/2), ...
                           [F_stop(i) F_stop2(i)]/(Fs/2), Rp, Rs);
        [b2{i}, a2{i}] = ellip(N, Rp, Rs, Wn, 'bandpass');
    end
    
    % Frequency response
    [H2{i}, w] = freqz(b2{i}, a2{i}, 1024, Fs);
    
    % Group delay
    Grp2{i} = grpdelay(b2{i}, a2{i}, 1024, Fs);
    
    % Impulse response
    h2{i} = impz(b2{i}, a2{i}, 120);
end

% Define colors for each filter
colors = lines(5); % Generates a set of distinct colors for plotting

% Combined Magnitude Spectrum Plot (Linear Scale)
figure(4);
hold on;
for i = 1:5
    plot(w, abs(H2{i}), 'Color', colors(i,:), 'DisplayName', ['Filter ', num2str(i)]);
end
hold off;
title('Combined Magnitude Response of Elliptic Filters (Linear Scale)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (Linear Scale)');
ylim([0, 1]); % Set y-axis range from 0 to 1
legend show;
grid on;

% Combined Phase Spectrum Plot with y-axis range -4 to 4
figure(5);
hold on;
for i = 1:5
    plot(w, angle(H2{i}), 'Color', colors(i,:), 'DisplayName', ['Filter ', num2str(i)]);
end
hold off;
title('Combined Phase Response of Elliptic Filters');
xlabel('Frequency (Hz)');
ylabel('Phase (Radians)');
ylim([-4, 4]); % Set y-axis range from -4 to 4
legend show;
grid on;

% Combined Group Delay Plot for Elliptic Filters
figure(6);
hold on;
for i = 1:5
    plot(w, Grp2{i}, 'Color', colors(i,:), 'DisplayName', ['Filter ', num2str(i)]);
end
hold off;
title('Combined Group Delay of Elliptic Filters');
xlabel('Frequency (Hz)');
ylabel('Group Delay (samples)');
legend show;
grid on;

% Plotting Impulse Response for each filter
figure(7);
for i = 1:5
    subplot(5, 1, i);
    stem(0:119, h2{i});
    title(['Filter ', num2str(i), ' Impulse Response']);
    xlabel('Samples'); ylabel('Amplitude');
    grid on;
end

% Plotting Pole-Zero plot for each filter
figure(8);
for i = 1:5
    subplot(3, 2, i);
    zplane(b2{i}, a2{i});
    title(['Filter ', num2str(i), ' Pole-Zero Plot']);
    grid on;
end



% Define gains for each filter
gains = [1, 2, 2.5, 5, 10];

% Initialize matrix to store filtered outputs and spectrograms
y = zeros(length(s), 5); % Filtered audio outputs
S = zeros(257, 1032, 5); % Spectrograms (257 frequency bins, 1032 time bins, 5 filters)

% Apply each filter with the specified gain
for i = 1:5
    % Filter the audio signal with the Chebyshev Type 1 filter
    filtered_audio = filter(b1{i}, a1{i}, s);
    
    % Apply the gain
    y(:, i) = gains(i) * filtered_audio;
    
    % Compute the spectrogram
    [S(:,:,i), f, t] = spectrogram(y(:, i), hanning(256), 128, 512, Fs);
end

% Compute and plot the spectrogram of the original audio
figure;
subplot(3, 2, 1);
spectrogram(s, hanning(256), 128, 512, Fs, 'yaxis');
title('Original Audio Spectrogram');
colorbar;

% Plot the spectrograms of the filtered audio signals
for i = 1:5
    subplot(3, 2, i+1);
    spectrogram(y(:, i), hanning(256), 128, 512, Fs, 'yaxis');
    title(['Filtered Audio Spectrogram - Filter ', num2str(i)]);
    colorbar;
end

% Optional: Save the spectrogram data for further analysis
% save('FilteredAudioData.mat', 'y', 'S');






end