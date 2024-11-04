function [result]=RMoGL(X,W,lambda1,lambda2,order,gt)

epson = 1e-7;
Clus_num = size(unique(gt),1); % number of clusters
V = size(X,2); % number of views
N = size(X{1},2);% number of data points
%normalized X
for i=1:V
    X{i} = NormalizeFea(X{i},0);  %Normalize feature
end
%Initilize Z,E,tensor G,multiplier Y,W
for i = 1:V
%     Z{i} = constructW_PKN(X{i},10);
    Z{i} = zeros(N,N);
    G{i} = zeros(N,N);
    B{i} = rand(size(X{i},1),size(W{i},1));
    W1{i} = zeros(N,N);  %multiplier
    W2{i} = zeros(N,N);  %multiplier
    E1{i} = zeros(size(X{i},1),N);
    E2{i} = zeros(N,N);
    C{i} = zeros(N,N);
    U = zeros(N,N);
    Y{i} = zeros(size(X{i},1),N);
end

weight = 1/V*ones(1,V);  


mu = 0.1;

mu_max = 10e10;
eta = 2;


Isconverg = 0;
iter = 1;

for v = 1:V
    Sn{v}=high_order_graph((X{v}+B{v}*W{v})',order);
end

%%  iteration
while(Isconverg == 0)
   
    %% update U
    theta_zv_sum = zeros(N,N);
    theta_sum = zeros(N,N);
    for i=1:V
        sn_sum = zeros(N,N);
        for j = 1:order
            sn_sum = sn_sum + Sn{i}{j};
        end
        theta_zv_sum = theta_zv_sum + weight(1,i)*(C{i}+sn_sum);
        theta_sum = theta_sum +  weight(1,i);
    end
    U_t = theta_zv_sum./theta_sum;

    Z0 = zeros(size(U_t));
    for is = 1:size(U_t,1)
        ind_c = 1:size(U_t,1);
        ind_c(is) = [];
        Z0(is,ind_c) = EProjSimplex_new(U_t(is,ind_c));
    end
    U = Z0;
    U = U - diag(diag(U));


    % == update E{i} ==
    for i=1:V
        F1 = X{i}+B{i}*W{i}-(X{i}+B{i}*W{i})*Z{i}+Y{i}./mu;
        F2 = Z{i}-C{i}+W2{i}./mu;
        F3 = [F1;F2];
        [Econcat] = solve_l1l2(F3,lambda2/mu);

        E1{i} = Econcat(1:size(F1,1),:);
        E2{i} = Econcat(size(F1,1)+1:size(F1,1)+size(F2,1),:);
    end


    for i=1:V
        [U1,S1,V1] = svd(Z{i}+W2{i}./mu);
        a1 = diag(S1)-lambda1/mu;
        a1(a1<0)=0;
        T1 = diag(a1);
        G{i} = U1*T1*V1';
    end



    % == update C{i} ==
    for i = 1:V
        sn_sum = zeros(N,N);
        for j = 1:order
            sn_sum = sn_sum + Sn{i}{j};
        end
        A1 = U - sn_sum;
        A2 = Z{i}-E2{i}+W1{i}./mu;
        C{i} = (2*weight(1,i)*A1+mu*A2)./(2*weight(1,i)+mu);
    end

    % == update B{i} ==
    for i=1:V
        C_1 = X{i}*Z{i}-X{i}+E1{i}-Y{i}./mu;
        D_1 = W{i}-W{i}*Z{i};
        B{i}= mu*C_1*D_1'*(inv(mu*D_1*D_1'));
    end

    % == update Z{i} ==
    for i = 1:V
        H{i} = X{i}+B{i}*W{i};
        A3 = H{i}-E1{i}+Y{i}./mu;
        A4 = C{i}+E2{i}-W1{i}./mu;
        A5 = G{i}-W2{i}./mu;
        Z{i} = inv(2*eye(N,N)+H{i}'*H{i})*(H{i}'*A3+A4+A5);
    end


    for v = 1:V
        Sn{v}=high_order_graph((X{v}+B{v}*W{v})',order);
    end


    % update weight
    for i = 1:V
        sn_sum = zeros(N,N);
        for j = 1:order
            sn_sum = sn_sum + Sn{i}{j};
        end
        weight(1,i) = 0.5/norm(U-(C{i}+sn_sum),'fro');
    end

    % == update Y{i} ==
    for i=1:V
        H2 = X{i}+B{i}*W{i};
        Y{i} = Y{i}+mu*(H2-H2*Z{i}-E1{i});
        W1{i} = W1{i}+mu*(Z{i}-C{i}-E2{i});
        W2{i} = W2{i}+mu*(Z{i}-G{i});
    end

    max_Z=0;
    max_Z_G=0;
    max_Z_C=0;
    Isconverg = 1;
    for k = 1:V
        H3 = X{k}+B{k}*W{k};
        if (norm(H3 - H3 * Z{k} - E1{k}, inf) > epson)
            history.norm_Z = norm(H3 - H3 * Z{k} - E1{k}, inf);
            Isconverg = 0;
            max_Z = max(max_Z, history.norm_Z);

        end

        if (norm(Z{k} - G{k}, inf) > epson)
            history.norm_Z_G = norm(Z{k} - G{k}, inf);
            Isconverg = 0;
            max_Z_G = max(max_Z_G, history.norm_Z_G);
        end

         if (norm(Z{k} - C{k}-E2{k}, inf) > epson)
            history.norm_Z_C = norm(Z{k} - C{k}-E2{k}, inf);
            Isconverg = 0;
            max_Z_C = max(max_Z_C, history.norm_Z_C);
        end

    end

    % == update  mu ==
    mu  = min(mu_max,mu*eta);
    iter = iter + 1;
    if (iter==50)
        Isconverg = 1;
    end

end


Clus = SpectralClustering((abs(U)+abs(U')/2),Clus_num); %
result1 = EvaluationMetrics(gt,Clus);
result = [result1(2),result1(3),result1(7)];
