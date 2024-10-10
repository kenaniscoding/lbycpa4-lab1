function [b1,b2,b3,b4,a1,a2,a3,a4,x,y,t,h1,h2,h3,h4,g1,g2,g3,g4,Same,s1,s2,s,so,Rejected,I1,I2,I3]=fname2
% De La Salle University
% Electronics and Computer Engineering Department
%
% Course        : LBYCPA4
% SECTION       : EQ1
% Submitted by  : Kenan Banal
% Submitted to  : Dr. Edwin Sybingco
%
% Exercise 2    : Discrete-time Systems
%
% Note: Check the instructions given in canvas

%% Task 1
%Encode the corresponding numerator coefficients and denominator 
%coefficients of the rational transfer functions of the causal LTI systems
% Use the following variables for the numerator coefficients: b1, b2, b3, and b4
% Use the following variables for the denominator coefficients: a1, a2, a3, and a4
% a.
b1 = [0.16, -0.48, 0.48, -0.16];
a1 = [1, 0.13, 0.52, 0.3];
% b.
b2 = [0.634, 0, -0.634];
a2 = [1, 0, -0.268];
% c.
b3 = [0.634, 0, 0.634];
a3 = [1, 0, 0.268];
% d.
b4 = [1, -5, 10];
a4 = [10, -5, 1];


%% Task 2
% Place your code here to generate the signal x
n=0:24;
wo =[0 pi/4 pi/2 3*pi/4 pi];

% Place your code here to generate the output sequences y
x = zeros(length(n), 5); 
for i = 1:5
    x(:, i) = cos(wo(i) * n);  
end
y = zeros(length(n), 5, 4);
y(:, :, 1) = filter(b1, a1, x);
y(:, :, 2) = filter(b2, a2, x);
y(:, :, 3) = filter(b3, a3, x);
y(:, :, 4) = filter(b4, a4, x);

% Place your code here to create the visualization of figure 1

figure(1);
for i = 1:5
    for j = 1:4
        subplot(5, 4, (i-1)*4 + j); % position 
        stem(n, x(:, i), 'b--', 'DisplayName', 'Input');  % Plot input
        hold on;
        stem(n, y(:, i, j), 'r', 'DisplayName', 'Output');  % Plot output
        hold off;
        title(['Sequence ', num2str(i), ' Filter ', num2str(j)]);
        xlabel('n');
        ylabel('Amplitude');
        legend('show');
    end
end
% Complete the code below to fillout the table
t=table([1;2;3;4],'VariableName',{'FilterNumber'});
t.FilterType(1)={'Highpass Filter'};
t.FilterType(2)={'Bandpass Filter'};
t.FilterType(3)={'Bandstop Filter'};
t.FilterType(4)={'Allpass Filter'};

%% Task 3: Impulse Response Plots
n = 0:24; % Time vector for 0 ≤ n ≤ 24
impulse = [1 zeros(1, 24)]; % Define impulse input
% Compute impulse response using filter for 0 <= n <= 24
h1 = filter(b1, a1, impulse);
h2 = filter(b2, a2, impulse);
h3 = filter(b3, a3, impulse);
h4 = filter(b4, a4, impulse);

% Using residuez for partial fraction expansion and impulse response
% Compute impulse responses using partial fraction expansion (residuez)
[r1, p1, k1] = residuez(b1, a1);
g1 = impz(r1, p1, k1, 25);

[r2, p2, k2] = residuez(b2, a2);
g2 = impz(r2, p2, k2, 25);

[r3, p3, k3] = residuez(b3, a3);
g3 = impz(r3, p3, k3, 25);

[r4, p4, k4] = residuez(b4, a4);
g4 = impz(r4, p4, k4, 25);

g1 = impz(b1, a1, n);
g2 = impz(b2, a2, n);
g3 = impz(b3, a3, n);
g4 = impz(b4, a4, n);

% Compare impulse responses h[n] and g[n]
Same = [isequal(h1, g1), isequal(h2, g2), isequal(h3, g3), isequal(h4, g4)];

% Plot the impulse responses
figure(2);
subplot(4, 2, 1); 
stem(n, h1, 'r'); 
title('h1[n] (Impulse Response)'); 
xlabel('n'); 
ylabel('Amplitude');

subplot(4, 2, 2); 
stem(n, g1, 'b'); 
title('g1[n] (Partial Fraction Expansion)'); 
xlabel('n'); 
ylabel('Amplitude');

subplot(4, 2, 3); 
stem(n, h2, 'r'); 
title('h2[n] (Impulse Response)'); 
xlabel('n'); 
ylabel('Amplitude');

subplot(4, 2, 4); 
stem(n, g2, 'b'); 
title('g2[n] (Partial Fraction Expansion)'); 
xlabel('n'); 
ylabel('Amplitude');

subplot(4, 2, 5); 
stem(n, h3, 'r'); 
title('h3[n] (Impulse Response)'); 
xlabel('n'); 
ylabel('Amplitude');

subplot(4, 2, 6); 
stem(n, g3, 'b'); 
title('g3[n] (Partial Fraction Expansion)'); 
xlabel('n'); 
ylabel('Amplitude');

subplot(4, 2, 7); 
stem(n, h4, 'r'); 
title('h4[n] (Impulse Response)'); 
xlabel('n'); 
ylabel('Amplitude');

subplot(4, 2, 8); 
stem(n, g4, 'b'); 
title('g4[n] (Partial Fraction Expansion)'); 
xlabel('n'); 
ylabel('Amplitude');
%% Task 4
n1 = 0:149;
% Generate s1, s2, and s
s1 = cos(n1 * pi / 10);
s2 = cos(n1 * pi / 4);
s = s1 + s2;

% Define the filter coefficients of the Notch filter
b = 0.9543 * [1, -2 * cos(pi / 4), 1];
a = [1, -1.9 * cos(pi / 4), 0.9025];

% Determine the output of the notch filter and store it in so
so = filter(b, a, s);  % Filtered output so[n1]

% Answer the question: Which of the sequences is rejected by the Notch filter (s1 or s2)?
Rejected = 's2 is rejected by the notch filter';
disp(Rejected);

% Create the visualization with stem plots
figure(3);

% Stem plot for s1[n1] and so[n1]
subplot(3, 1, 1);
stem(n1, s1, 'b--', 'DisplayName', 's_1[n_1]');
hold on;
stem(n1, so, 'r', 'DisplayName', 's_o[n_1]');
hold off;
title('s_1[n_1] overlayed with s_o[n_1]');
xlabel('n_1');
ylabel('Amplitude');
legend('show');

% Stem plot for s2[n1] and so[n1]
subplot(3, 1, 2);
stem(n1, s2, 'b--', 'DisplayName', 's_2[n_1]');
hold on;
stem(n1, so, 'r', 'DisplayName', 's_o[n_1]');
hold off;
title('s_2[n_1] overlayed with s_o[n_1]');
xlabel('n_1');
ylabel('Amplitude');
legend('show');

% Stem plot for s[n1] and so[n1]
subplot(3, 1, 3);
stem(n1, s, 'b--', 'DisplayName', 's[n_1]');
hold on;
stem(n1, so, 'r', 'DisplayName', 's_o[n_1]');
hold off;
title('s[n_1] overlayed with s_o[n_1]');
xlabel('n_1');
ylabel('Amplitude');
legend('show');

%% Task 5
% Retrieve lena_gray.tiff
I = imread("lena_gray.tiff");

% Place you code here to represent the three filters in terms of their
% convolution kernels. You can use any variables to represents the three
% filters.
h(:,:,1) = [-1 -1 -1; -1 8 -1; -1 -1 -1];
h(:,:,2) = [-1 -1 -1; -1 9 -1; -1 -1 -1];
h(:,:,3) = [0 1 2; -1 0 1; -2 -1 0];
% Place your code here to determine I1. I1 should be uint8
I1 = uint8(filter2(h(:,:,1), I, 'same'));
% Place your code here to determine I2. I2 should be uint8
I2 = uint8(filter2(h(:,:,2), I, 'same'));
% Place your code here to determine I3. I3 should be uint8
I3 = uint8(filter2(h(:,:,3), I, 'same'));

% Place your code here to create the visualization of figure 4.
figure(4)
subplot(2, 2, 1);
imshow(I);
title('Original Image');
subplot(2, 2, 2);
imshow(I1);
title('Filtered Image I1');
subplot(2, 2, 3);
imshow(I2);
title('Filtered Image I2');
subplot(2, 2, 4);
imshow(I3);
title('Filtered Image I3');
end
