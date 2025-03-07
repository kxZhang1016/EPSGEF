function S = f_L1(X, h_input_2, v_input_2, lambda, gamma)

% ȫ�ֳ�����Ĭ��ֵ
% image smoothing
MAX_ITER = 100; % �������� image smoothing
rho = 2; % �������ճͷ����� 3.2 1.5

% texture smoothing
% MAX_ITER = 200; % texture smoothing
% rho = 2; % �������ճͷ����� 2.3

% eliminate noise
% MAX_ITER = 100; % 40
% rho = 1.8; % 2.5

% ����Ԥ����
X = im2double(X);
[row, col, cha] = size(X);
w = gamma / rho; % w����shrink����

% ADMM Solver
S = zeros(row, col, cha);
% S = X;
Lx = zeros(row, col, cha);
Ly = zeros(row, col, cha);
dx = zeros(row, col, cha);
dy = zeros(row, col, cha);

fx = [1, -1]; % ˮƽ���
fy = [1; -1]; % ��ֱ���
sizeI2D = [row, col];
otfFx = psf2otf(fx, sizeI2D); % ��ˮƽ�����ϵ�һά����ɢ����(PSF)ת��Ϊ��Ӧ�Ķ�ά��ѧ���ݺ���(OTF)
otfFy = psf2otf(fy, sizeI2D);
Normin = fft2(X);
Denormin1 = abs(otfFx).^2 + abs(otfFy).^2;
if cha > 1
    Denormin1 = repmat(Denormin1, [1,1,cha]);
end
Denormin = 1 + (lambda + rho/2) * Denormin1;

for k = 1:MAX_ITER % ����
    
    %S-update
    h1 = Lx - dx;
    v1 = Ly - dy;
    h2 = h_input_2;
    v2 = v_input_2;
    Normin11 = [h1(:,end,:) - h1(:, 1,:), -diff(h1,1,2)];
    Normin12 = Normin11 + [v1(end,:,:) - v1(1, :,:); -diff(v1,1,1)]; % g*
    Normin21 = [h2(:,end,:) - h2(:, 1,:), -diff(h2,1,2)];
    Normin22 = Normin21 + [v2(end,:,:) - v2(1, :,:); -diff(v2,1,1)];
    FS = (Normin + lambda * fft2(Normin22) + (rho/2) * fft2(Normin12)) ./ Denormin; % FFT
    S = real(ifft2(FS));  % IFFT
    
    % L-update �����������ķ�ʽ����
    h = [diff(S,1,2), S(:,1,:) - S(:,end,:)];
    v = [diff(S,1,1); S(1,:,:) - S(end,:,:)];
    Ax_hat1 = h + dx;
    Ax_hat2 = v + dy;
    Lx = max(abs(Ax_hat1) - w, 0).*sign(Ax_hat1);
    Ly = max(abs(Ax_hat2) - w, 0).*sign(Ax_hat2);
    
    % d-update
    dx = (h - Lx) * rho + dx;
    dy = (v - Ly) * rho + dy;
    
end

end

