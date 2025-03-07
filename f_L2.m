function S = f_L2(X,h_input_2,v_input_2,lambda,gamma)

X = im2double(X); 
[row, col, cha] = size(X); 

fx = [1, -1];
fy = [1; -1];
sizeI2D = [row,col];
otfFx = psf2otf(fx,sizeI2D);
otfFy = psf2otf(fy,sizeI2D);
Normin = fft2(X); % FFT
Denormin1 = abs(otfFx) .^ 2 + abs(otfFy ) .^ 2; 
if cha > 1 
    Denormin1 = repmat(Denormin1, [1,1,cha]); 
end
Denormin = 1 +  (lambda + gamma) * Denormin1; 

h = h_input_2;
v = v_input_2;
Normin1 = [h(:,end,:) - h(:, 1,:), -diff(h,1,2)]; 
Normin2 = Normin1 + [v(end,:,:) - v(1, :,:); -diff(v,1,1)]; 
FS = (Normin + lambda * fft2(Normin2)) ./ Denormin; 
S = real(ifft2(FS));  

end

