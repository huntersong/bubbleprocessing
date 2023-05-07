function [xx_selected,yy_selected] = bezier_points(P0,P1,P2) 
    
    t = linspace(0, 1, 1000);
    x = (1-t).^2 * P0(1,1) + 2*(1-t).*t*P1(1,1) + t.^2*P2(1,1);
    y = (1-t).^2 * P0(1,3) + 2*(1-t).*t*P1(1,3) + t.^2*P2(1,3);
    
    idx = find(x>=P0(1,1) & x<=P2(1,1));% 提取x=25到x=100之间的点
    x_selected = round(x(idx)); % 四舍五入为整数
    y_selected = round(y(idx));
    xx_selected(1,1) = x_selected(1,1);  yy_selected(1,1) = y_selected(1,1);
    m=1;
    for n = 1:(length(x_selected)-1)
        if x_selected(1,n) ~= x_selected(1,n+1)
            m = m +1;
            xx_selected(1,m) = x_selected(1,n+1);
            yy_selected(1,m) = y_selected(1,n+1);
        end
    end
%plot(xx_selected, yy_selected, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
%axis equal;
end
