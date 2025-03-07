clear;

hdr = double(hdrread('7.hdr'));
figure;imshow(hdr);
I = 0.2989*hdr(:,:,1) + 0.587*hdr(:,:,2) + 0.114*hdr(:,:,3);
logI = log(I+eps);
X = logI;

lambda = 10; 
gamma = 10;  
sigma = 0.12;  
alpha = 2;  
sigma_s = 16;  
sigma_r = 0.01;  

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

base = Em_BLF_L1_ColdStart(X, h_input_2, v_input_2, lambda, gamma);

compression = 0.25; 
detail = X - base;
OUT = base * compression +  1.2 * detail;
OUT = exp(OUT);

OUT = OUT./ I;
OUT = hdr .* padarray(OUT, [0 0 2], 'circular' , 'post');

gamma = 1.0/2.2;
bias = -min(OUT(:));
gain = 0.45;  
OUT_SP = (gain * (OUT + bias)).^gamma;

figure;
imshow(OUT_SP);




