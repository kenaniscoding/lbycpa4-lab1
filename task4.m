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

