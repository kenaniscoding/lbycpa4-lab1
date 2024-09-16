function [x,n,NRate,s1,s2,s3,yq,MSE,Response,I]=fname1

% De La Salle University
% Electronics and Computer Engineering Department
%
% Course        : LBYCPA4
% SECTION       : EQ1
% Submitted by  : Kenan Banal
% Submitted to  : Dr. Edwin Sybingco
%
% Exercise 1: Discrete-time Signals
%
% Note: Check the instructions given in the canvas

%% Task 1
% Place your code here to generate x1, x2, x3,x4,x5, and x6
n=-25:25;
delta = @(n) double(n==0);
u = @(n) double(n>=0);
x1 = -delta(n+3)+delta(n)+delta(n-5);
x2 = u(n+3)-u(n-5);
x3 = ((-27/40).^n).*u(n);
x4 = cos(((3*pi/12)*n)-(pi/6));
x5 = ((-27/40).^n).*cos(((3*pi/12)*n)-(pi/6)).*u(n);
x6 = cos(((3*pi/12)*n)-(pi/6)).*(u(n+3)-u(n-5));

% Place your code here to create the visualization of Figure 1 
figure;
subplot(3, 2, 1);
stem(n, x1, 'filled');
title('x1[n]');
xlabel('n');
ylabel('x1[n]');
 
subplot(3, 2, 2);
stem(n, x2, 'filled');
title('x2[n]');
xlabel('n');
ylabel('x2[n]');
 
subplot(3, 2, 3);
stem(n, x3, 'filled');
title('x3[n]');
xlabel('n');
ylabel('x3[n]');
 
subplot(3, 2, 4);
stem(n, x4, 'filled');
title('x4[n]');
xlabel('n');
ylabel('x4[n]');
 
subplot(3, 2, 5);
stem(n, x5, 'filled');
title('x5[n]');
xlabel('n');
ylabel('x5[n]');
 
subplot(3, 2, 6);
stem(n, x6, 'filled');
title('x6[n]');
xlabel('n');
ylabel('x6[n]');
% The following code will concatenate the six discrete-time sequences. x
% represents the disctrete-time sequences with 6 channels.
x=[x1' x2' x3' x4' x5' x6'];
%% Task 2
% Define the time variable and the original signal
t = -0.5:0.001:0.5;
st = cos(40*pi*t).*cos(2*pi*t);
NRate = 42;  % Nyquist rate in Hz

% Define Fs1, t1, and s1
Fs1 = 5*NRate;
t1 = -0.5:1/Fs1:0.5;
s1 = cos(40*pi*t1).*cos(2*pi*t1);

% Define Fs2, t2, and s2
Fs2 = NRate;
t2 = -0.5:1/Fs2:0.5;
s2 = cos(40*pi*t2).*cos(2*pi*t2);

% Define Fs3, t3, and s3
Fs3 = (3/4)*NRate;
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
    
    sound(2^-(nbits(i)-1)*yq(:,i), Fs);
    pause(2);
end
% Compute the mean square error
MSE = zeros(1, length(nbits));  
for i = 1:length(nbits)
   MSE(i) = sum((yq(:,i) - y).^2) / N;
end

% plot the mean square error in figure 4. Label the corresponding graph.
figure(4);
plot(nbits, MSE, '-o', 'LineWidth', 2);
xlabel('Number of Bits');
ylabel('Mean Square Error');
grid on;
% Is the MSE consistent with the played back audio? (Yes or No). Store your
% response as a string in variable Response.
Response = 'Yes';
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
