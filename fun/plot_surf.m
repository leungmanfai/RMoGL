function [] = plot_surf(str_x,str_y,y_nama_,x_name,z,acc_nmi,name)
figure('Position', [200, 200, 820, 400]);
Go = surf(z);
%colormap('winter');
colormap('cool');
%beta:
set(gca,'Yticklabel', str_x);
set(gca, 'YTick', 1:10);
ylabel(y_nama_,'FontWeight', 'bold');

%gamma:
set(gca,'Xticklabel', str_y);
set(gca, 'XTick', 1:10);
xlabel(x_name,'FontWeight', 'bold');

zlim([0,1]);


c=colorbar;
c.Location='NorthOutside';
%c.Position(1) = 0.042;  % 调整左边距
%c.Position(3) = 0.011;    % 调整宽度，设置为您需要的值

ax = gca;
ax.FontSize = 24;  % 设置为您需要的字体大小
ax.XLabel.FontSize = 36;  % 设置为您需要的字体大小
ax.YLabel.FontSize = 36;  % 设置为您需要的字体大小
% 设置坐标轴刻度标签的字体大小




if acc_nmi=="acc"
    zlabel('ACC');
    ax = gca;  % 获取当前的坐标轴句柄
    ax.LineWidth = 1.2;  % 设置坐标轴线的宽度
    print(name, '-depsc');
end
if acc_nmi=="nmi"
    zlabel('NMI');
    ax = gca;  % 获取当前的坐标轴句柄
    ax.LineWidth = 1.2;  % 设置坐标轴线的宽度
    print(name, '-depsc');
end

end