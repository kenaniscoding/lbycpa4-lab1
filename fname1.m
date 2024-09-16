function [x,n,NRate,s1,s2,s3,yq,MSE,Response,I]=fname1

% De La Salle University
% Electronics and Computer Engineering Department
%
% Course        : LBYCPA4
% SECTION       :
% Submitted by  :
% Submitted to  : Dr. Edwin Sybingco
%
% Exercise 1: Discrete-time Signals
%
% Note: Check the instructions given in the canvas

%% Task 1
% Place your code here to generate x1, x2, x3,x4,x5, and x6
n=-25:25;

% Place your code here to create the visualization of Figure 1 

% The following code will concatenate the six discrete-time sequences. x
% represents the disctrete-time sequences with 6 channels.
x=[x1' x2' x3' x4' x5' x6'];
%% Task 2
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
title(['Sampled Rate at Fs1 = ', num2str(Fs1), ' Hz']);
xlabel('Time (s)');
ylabel('Amplitude');
grid on; hold off;

% Plot sampled signal s2 with Fs2
subplot(3,1,2);
plot(t, st, 'r'); hold on;
stem(t2, s2, 'b');
title(['Sampled Rate at Fs2 = ', num2str(Fs2), ' Hz']);
xlabel('Time (s)');
ylabel('Amplitude');
grid on; hold off;

% Plot sampled signal s3 with Fs3
subplot(3,1,3);
plot(t, st, 'r'); hold on;
stem(t3, s3, 'b');
title(['Sampled Rate at Fs3 = ', num2str(Fs3), ' Hz']);
xlabel('Time (s)');
ylabel('Amplitude');
grid on; hold off;

% Adjust the layout for better readability
sgtitle('Task 2');

%% Task 3: Quantization of audio into n number of bits

% read the audio tononoise8mono.wav and store audio data and sampling rate
% in y and Fs, respectively

% Store the number of samples to N and generate the corresponding
% continuous time t

% Convert the audio data into 8-bit signed integer


% Quatized the audio data into n-bits. Store the quantized audio as ynbits.
% Example y7 for 7 bits, y6 for 6 bits and so on.Store the audio data in yq 
% with the following components: yq=[y8 y7...y1]
nbits=7:-1:1;


% Plot the audio data in figure 3. Label the graph

% Compute the mean square error

% plot the mean square error in figure 4. Label the corresponding graph.

% Is the MSE consistent with the played back audio? (Yes or No). Store your
% response as a string in variable Response.

%% Task 4: Image Generation
% Place your code here to create the image I
% Initialize I
I = zeros(140,140,3);
bw = 20; %bar width

% Generate white color
I(:, 1:bw, :) = repmat(reshape([255, 255, 255], 1, 1, 3), 140, bw);

% Generate Yellow color 
I(:, bw+1:2*bw, :) = repmat(reshape([255, 255, 0], 1, 1, 3), 140, bw);

% Generate Cyan color 
I(:, 2*bw+1:3*bw, :) = repmat(reshape([0, 255, 255], 1, 1, 3), 140, bw);

% Generate Green color 
I(:, 3*bw+1:4*bw, :) = repmat(reshape([0, 255, 0], 1, 1, 3), 140, bw);

% Generate Magenta color 
I(:, 4*bw+1:5*bw, :) = repmat(reshape([255, 0, 255], 1, 1, 3), 140, bw);

% Generate Red color 
I(:, 5*bw+1:6*bw, :) = repmat(reshape([255, 0, 0], 1, 1, 3), 140, bw);

% Generate Blue color
I(:, 6*bw+1:end, :) = repmat(reshape([0, 0, 255], 1, 1, 3), 140, bw);

% Display the generated image in figure 5
figure(5);
imshow(I);
title('Task 4');


end
