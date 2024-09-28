function F = W2F(W)

% 将代码X+EW中的W（1代表缺失）转化为F，F{i}为Nx1，0代表缺失。

% load("ORL_p_0.1_balaced_incomplete.mat");
% clear all;
% W{1} = [1 0 0 0 0;
%      0 0 1 0 0
%      0 1 0 0 0];

% 示例矩阵 A，你可以替换成你的实际矩阵
N = size(W{1},2);
V = size(W,2);
for i = 1:V
    A = W{i};
    % 找到矩阵 A 中值为 1 的元素的列索引
    [col_indices, ~] = find(A' == 1);

    % 创建新矩阵 B
    B = col_indices';
    C = ones(N,1);
    C(B) = 0;
    F(:,i) = C;
end
end
