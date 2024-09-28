function W2 = high_order_graph(x, order)
% 计算余弦相似性矩阵
N = size(x,1);
if order ==1
    S = constructW_PKN(x',10);
    W2{1} = S;
    W2{1} = W2{1} - diag(diag(W2{1}));
else
   % S = 1 - pdist2(x, x, 'cosine');
    S = constructW_PKN(x',10);
    W2{1} = S;
    W2{1} = W2{1} - diag(diag(W2{1}));

    for o = 2:(order)
        W2{o} = W2{o-1}*W2{1};
        W2{o} = W2{o} - diag(diag(W2{o}));
    end
end


