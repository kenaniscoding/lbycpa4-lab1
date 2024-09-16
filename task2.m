%% Task 2: Nyquist Rate and Signal Visualization with Subplots

% Define the time variable and the original signal
t = -0.5:0.001:0.5;
st = cos(40*pi*t).*cos(2*pi*t);
fc = 42;  % Nyquist rate in Hz
%(1/2)*(cos(38*pi*t2)) + (1/2)*(cos(42*pi*t2));

% Define Fs1, t1, and s1
Fs1 = 5*fc;
t1 = -0.5:1/Fs1:0.5;
s1 = cos(40*pi*t1).*cos(2*pi*t1);

% Define Fs2, t2, and s2
Fs2 = fc;
t2 = -0.5:1/Fs2:0.5;
s2 = cos(40*pi*t2).*cos(2*pi*t2);

% Define Fs3, t3, and s3
Fs3 = (3/4)*fc;
t3 = -0.5:1/Fs3:0.5;
s3 = cos(40*pi*t3).*cos(2*pi*t3);

% Create subplots for visualization
figure(2);

% Plot sampled signal s1 with Fs1
subplot(3,1,1);
plot(t, st, 'r'); hold on;
stem(t1, s1, 'b');
title(['Sampled Rate at Fs1 = ', num2str(Fs1)]);
xlabel('Time (s)');
ylabel('Amplitude');
grid on; hold off;

% Plot sampled signal s2 with Fs2
subplot(3,1,2);
plot(t, st, 'r'); hold on;
stem(t2, s2, 'b');
title(['Sampled Rate at Fs2 = ', num2str(Fs2)]);
xlabel('Time (s)');
ylabel('Amplitude');
grid on; hold off;

% Plot sampled signal s3 with Fs3
subplot(3,1,3);
plot(t, st, 'r'); hold on;
stem(t3, s3, 'b');
title(['Sampled Rate at Fs3 = ', num2str(Fs3)]);
xlabel('Time (s)');
ylabel('Amplitude');
grid on; hold off;

% Adjust the layout for better readability
sgtitle('Task 2');
