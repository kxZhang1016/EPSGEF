function S = f_Lp(X, h_input_2, v_input_2, lambda, gamma)

MAX_ITER = 5; 
rho = 2; 
p = 0.7; 
smallNum = 0.0001;

X = im2double(X);
[row, col] = size(X); 
a = row * col;
w = gamma / rho; 

% ADMM Solver
S = zeros(row, col); 
Lx = zeros(row, col); 
Ly = zeros(row, col);
dx = zeros(row, col);
dy = zeros(row, col);
weight_x = ones(row, col); 
weight_y = ones(row, col);
I = ones(row, col);

for k = 1:MAX_ITER
    weight_x_sparse = spdiags(weight_x(:), 0, a, a);
    weight_y_sparse = spdiags(weight_y(:), 0, a, a);
    
    h1 = reshape(weight_x_sparse' * (Lx(:) - dx(:)), row, col);
    v1 = reshape(weight_y_sparse' * (Ly(:) - dy(:)), row, col);
    h2 = h_input_2;
    v2 = v_input_2;
    
    Normin11 = [h1(:,end) - h1(:, 1), -diff(h1,1,2)];
    Normin12 = Normin11 + [v1(end,:) - v1(1, :); -diff(v1,1,1)];
    Normin21 = [h2(:,end) - h2(:, 1), -diff(h2,1,2)];
    Normin22 = Normin21 + [v2(end,:) - v2(1, :); -diff(v2,1,1)];

    A = X(:) + lambda * Normin22(:) + (rho/2) * Normin12(:);
    
    [Dx, Dy] = Diff(X);
    B = spdiags(I(:), 0, a, a) + lambda * (Dx' * Dx + Dy' * Dy) + (rho/2) * (Dx' * (weight_x_sparse' * weight_x_sparse) * Dx + Dy' * (weight_y_sparse' * weight_y_sparse) * Dy);
    
    if exist('ichol', 'builtin')
        L = ichol(B, struct('michol', 'on'));    
        tin = A;
        [tout, flag] = pcg(B, tin(:), 0.1, 100, L, L'); 
        S = reshape(tout, row, col);  
    else
        tin = A;
        tout = B \ tin(:);
        S = reshape(tout, row, col);  
    end
    
    h = [diff(S, 1, 2), S(:, 1) - S(:, end)];
    v = [diff(S, 1, 1); S(1, :) - S(end, :)];
    Ax_hat1 = weight_x .* h + dx;
    Ax_hat2 = weight_y .* v + dy;
    Lx = max(abs(Ax_hat1) - w, 0) .* sign(Ax_hat1);
    Ly = max(abs(Ax_hat2) - w, 0) .* sign(Ax_hat2);
    
    dx = (weight_x .* h - Lx) * rho + dx;
    dy = (weight_y .* v - Ly) * rho + dy;
    
    weight_x = abs([diff(S, 1, 2), S(:, 1) - S(:, end)] + smallNum) .^ (p - 1);
    weight_y = abs([diff(S, 1, 1); S(1, :) - S(end, :)] + smallNum) .^ (p - 1);
end    

end

function [Cx, Cy] = Diff(I)
    [r, c] = size(I);
    x = ones(r, c);
    k = r * c;
    x = x(:);
    y = -1 * x(:);
    B(:,1) = y;
    B(:,2) = x;
    B(:,3) = y;

    Cy = spdiags(B, [-k+1, 0, 1], k, k);
    
    D(:,1) = y;
    D(:,2) = x;
    D(:,3) = y;
    Cx = spdiags(D, [r-k, 0, r], k, k);
end
