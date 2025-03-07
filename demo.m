clear;

Img = im2double(imread('D:\DataSet\BSDS500\data\images\val\145086.jpg'));

X = Img;

% L0.7 image smoothing
% lambda = 3; 
% gamma = 0.05;  % smooth 0.2 denoise 0.5  
% sigma = 0.05;  
% alpha = 2;  
% sigma_s = 12;  
% sigma_r = 0.1;


% % L1 image smoothing
% lambda = 5; 
% gamma = 1;  
% sigma = 0.12;  
% alpha = 2;  
% sigma_s = 16;  
% sigma_r = 0.01;  

% L1 detail enhancement
% lambda = 5; 
% gamma = 0.3;  
% sigma = 0.12;  
% alpha = 2;  
% sigma_s = 16;  
% sigma_r = 0.01;  

% L1 detail enhancement
% lambda = 1;  
% gamma = 0.4;  
% sigma = 0.8;  
% alpha = 3;  
% sigma_s = 12; 
% sigma_r = 0.3;  


% L1 cartoon image artificial remove
% lambda = 0.5;  
% gamma = 0.1;  
% sigma = 1;  
% alpha = 2;  
% sigma_s = 16;
% sigma_r = 0.01;  

% L1 texture smoothing-test
% lambda = 0.1;  
% gamma = 0.8;  
% sigma = 0.4;  
% alpha = 3;  
% sigma_s = 16; 
% sigma_r = 0.01;  

% L1 texture smoothing
lambda = 0.1;  
gamma = 0.8;  
sigma = 0.4;  
alpha = 3;  
sigma_s = 16; 
sigma_r = 0.01;  

% L1 eliminate noise
% lambda = 0.6; 
% gamma = 1.8;  
% sigma = 0.6;  
% alpha = 5; 
% sigma_s = 14;  
% sigma_r = 0.02;  

% L2 detail enhancing
% lambda = 2.5;
% gamma = 0.3;
% sigma = 0.2;
% alpha = 2;
% sigma_s = 16;
% sigma_r = 0.001;

% L2 image smoothing
% lambda = 100; 
% gamma = 100; 
% sigma = 0.12; 
% alpha = 2; 
% sigma_s = 16; 
% sigma_r = 0.02; 

% L2 texture smoothing
% lambda = 100; 
% gamma = 100;  
% sigma = 0.4; 
% alpha = 5; 
% sigma_s = 14; 
% sigma_r = 0.01;

% L2 noise smoothing
% lambda = 200; 
% gamma = 80;  
% sigma = 0.6; 
% alpha = 5; 
% sigma_s = 16; 
% sigma_r = 0.02;


[row, col, cha] = size(X);

h_input_1 = [diff(X, 1, 2), X(:, 1, :) - X(:, end, :)];
v_input_1 = [diff(X, 1, 1); X(1, :, :) - X(end, :, :)];

tt1 = (abs(h_input_1)>=sigma);
tt2 = (abs(v_input_1)>=sigma);

hh_input = h_input_1(tt1);
hv_input = v_input_1(tt2);

h_input_1 = sign(h_input_1)*sigma.*(abs(h_input_1)/sigma).^alpha; 
v_input_1 = sign(v_input_1)*sigma.*(abs(v_input_1)/sigma).^alpha;

h_input_1(tt1) = hh_input;
v_input_1(tt2) = hv_input;
 
h_input_1(:, end, :) = X(:, 1, :) - X(:, end, :);
v_input_1(end, :, :) = X(1, :, :) - X(end, :, :);

h_input_2 = zeros(row, col, cha);
v_input_2 = zeros(row, col, cha);

if cha == 1
    h_input_2 = bilateralFilter(h_input_1, h_input_1, min(h_input_1(:)), max(h_input_1(:)), sigma_s, sigma_r);
    v_input_2 = bilateralFilter(v_input_1, v_input_1, min(v_input_1(:)), max(v_input_1(:)), sigma_s, sigma_r);
else
    for k = 1:cha
        I_h = h_input_1(:, :, k);
        h_input_2(:, :, k) = bilateralFilter(I_h, I_h, min(I_h(:)), max(I_h(:)), sigma_s, sigma_r);
    
        I_v = v_input_1(:, :, k);
        v_input_2(:, :, k) = bilateralFilter(I_v, I_v, min(I_v(:)), max(I_v(:)), sigma_s, sigma_r);
    end
end

% if cha == 1
%     S = f_Lp(X, h_input_2, v_input_2, lambda, gamma);
% else
%     for i = 1 : cha
%         S(:, :, i) = f_Lp(X(:, :, i), h_input_2(:, :, i), v_input_2(:, :, i), lambda, gamma);
%     end
% end

S = f_L1(X, h_input_2, v_input_2, lambda, gamma);

figure;
imshow(S);