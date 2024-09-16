%% Task 1
% Place your code here to generate x1, x2, x3, x4, x5, and x6
n=-25:25;

%x1
d1=myunit(n+3);
d2=myunit(n);
d3=myunit(n-5);
x1=-d1+d2+d3;

%x2
u1=mystep(n+3);
u2=mystep(n-5);
x2=u1-u2;

%x3
u3=mystep(n);
x3=((-27/40).^n).*u3;

%x4
x4=cos(((3*pi)/12)*n-pi/6);

%x5
x5=((-27/40).^n).*cos((3*pi)/12*n-pi/6).*u3;

%x6
x6=cos((3*pi)/12*n-pi/6).*(u1-u2);

% Visualization of Figure 1 
figure(1)

%x1
subplot(3,2,1);
stem(n,x1)
grid on
xlabel('n')
ylabel('x1[n]')
title('Sequence x1')

%x2
subplot(3,2,2);
stem(n,x2)
grid on
xlabel('n')
ylabel('x2[n]')
title('Sequence x2')

%x3
subplot(3,2,3);
stem(n,x3)
grid on
xlabel('n')
ylabel('x3[n]')
title('Sequence x3')

%x4
subplot(3,2,4);
stem(n,x4)
grid on
xlabel('n')
ylabel('x4[n]')
title('Sequence x4')

%x5
subplot(3,2,5);
stem(n,x5)
grid on
xlabel('n')
ylabel('x5[n]')
title('Sequence x5')

%x6
subplot(3,2,6);
stem(n,x6)
grid on
xlabel('n')
ylabel('x6[n]')
title('Sequence x6')

% The following code will concatenate the six discrete-time sequences. x
% represents the discrete-time sequences with 6 channels.
x=[x1' x2' x3' x4' x5' x6'];

% Custom function definitions
function y = myunit(n)
    y = zeros(1, length(n));
    [a, b] = find(n == 0);
    y(b) = 1;
end

function y = mystep(n)
    y = zeros(1, length(n));
    [a, b] = find(n >= 0);
    y(b) = 1;
end

