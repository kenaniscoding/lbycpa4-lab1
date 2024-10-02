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
        plot(n, x(:, i), 'b--', 'DisplayName', 'Input');  % Plot input
        hold on;
        plot(n, y(:, i, j), 'r', 'DisplayName', 'Output');  % Plot output
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

%% Task 3
n = 0:24; % Time vector for 0 ≤ n ≤ 24
impulse = [1 zeros(1, 24)]; % Define impulse input

% Define filter coefficients for four filters
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


h1 = filter(b1, a1, impulse);
h2 = filter(b2, a2, impulse);
h3 = filter(b3, a3, impulse);
h4 = filter(b4, a4, impulse);


[r1, p1, k1] = residuez(b1, a1);
%g1
g1 = zeros(1, 25);
for i = 1:length(r1)
    g1 = g1 + r1(i) * (p1(i) .^ n);
end
g1 = g1 + polyval(k1, n);

[r2, p2, k2] = residuez(b2, a2);
%g2
g2 = zeros(1, 25);
for i = 1:length(r2)
    g2 = g2 + r2(i) * (p2(i) .^ n);
end
g2 = g2 + polyval(k2, n);

[r3, p3, k3] = residuez(b3, a3);
%g3
g3 = zeros(1, 25);
for i = 1:length(r3)
    g3 = g3 + r3(i) * (p3(i) .^ n);
end
g3 = g3 + polyval(k3, n);

[r4, p4, k4] = residuez(b4, a4);
%g4
g4 = zeros(1, 25);
for i = 1:length(r4)
    g4 = g4 + r4(i) * (p4(i) .^ n);
end
g4 = g4 + polyval(k4, n);


Same = zeros(1, 4);
Same(1) = isequal(h1, g1);
Same(2) = isequal(h2, g2);
Same(3) = isequal(h3, g3);
Same(4) = isequal(h4, g4);


subplot(4, 2, 1); 
plot(n, h1, 'r', 'LineWidth', 1); 
title('h1[n]'); 
xlabel('n'); 
ylabel('Impulse Response');

subplot(4, 2, 2); 
plot(n, g1, 'b', 'LineWidth', 1); 
title('g1[n]'); 
xlabel('n'); 
ylabel('Impulse Response');

subplot(4, 2, 3); 
plot(n, h2, 'r', 'LineWidth', 1); 
title('h2[n]'); 
xlabel('n'); 
ylabel('Impulse Response');

subplot(4, 2, 4); 
plot(n, g2, 'b', 'LineWidth', 1); 
title('g2[n]'); 
xlabel('n'); 
ylabel('Impulse Response');

subplot(4, 2, 5); 
plot(n, h3, 'r', 'LineWidth', 1); 
title('h3[n]'); 
xlabel('n'); 
ylabel('Impulse Response');

subplot(4, 2, 6); 
plot(n, g3, 'b', 'LineWidth', 1); 
title('g3[n]'); 
xlabel('n'); 
ylabel('Impulse Response');

subplot(4, 2, 7); 
plot(n, h4, 'r', 'LineWidth', 1); 
title('h4[n]'); 
xlabel('n'); 
ylabel('Impulse Response');

subplot(4, 2, 8); 
plot(n, g4, 'b', 'LineWidth', 1); 
title('g4[n]'); 
xlabel('n'); 
ylabel('Impulse Response');

%% Task 4
n1 = 0:149;
% Place your code here to generate s1, s2, and s
s1=cos(n1*pi/10);
s2=cos(n1*pi/4);
s=s1+s2;
% place your code here define the filter coefficients of the Notch filter
b = 0.9543 * [1, -2*cos(pi/4), 1];
a = [1, -1.9*cos(pi/4), 0.9025];
% Place your code here to determine the output of the notch filter and
% store it in so
so = filter(b, a, s);  % Filtered output so[n1]
figure(3);
% Plot s1[n1] and so[n1]

% Answer the question: Which of the sequences is rejected by the Notch
% filter (s1 or s2)? Store your response as a string using the variable
% Rejected
Rejected = 's2 is rejected by the notch filter';
disp(Rejected);
% Place your code here to create the visualization of figure 3.
subplot(3, 1, 1);
plot(n1, s1, 'b--', 'DisplayName', 's_1[n_1]');
hold on;
plot(n1, so, 'r', 'DisplayName', 's_o[n_1]');
hold off;
title('s_1[n_1] overlayed with s_o[n_1]');
xlabel('n_1');
ylabel('Amplitude');
legend('show');
% Plot s2[n1] and so[n1]
subplot(3, 1, 2);
plot(n1, s2, 'b--', 'DisplayName', 's_2[n_1]');
hold on;
plot(n1, so, 'r', 'DisplayName', 's_o[n_1]');
hold off;
title('s_2[n_1] overlayed with s_o[n_1]');
xlabel('n_1');
ylabel('Amplitude');
legend('show');
% Plot s[n1] and so[n1]
subplot(3, 1, 3);
plot(n1, s, 'b--', 'DisplayName', 's[n_1]');
hold on;
plot(n1, so, 'r', 'DisplayName', 's_o[n_1]');
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
