clear all;

%rng(0);

load("YaleA_0.2_unbalanced_incomplete.mat");

lambda1 =0.01 ;
lambda2 = 0.1;
order = 2;
result = RMoGL(X,W,lambda1,lambda2,order,gt);
fprintf("NMI = %5.4f，Purity = %5.4f，ARI = %5.4f\n",result(1),result(2),result(3));


