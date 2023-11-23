
point1_img_camera1 = [25,147,-1000];  point2_img_camera1 = [238,147,-1100];  ...
        point1_img_camera2 = [100,147,-900];  point2_img_camera2 = [80,147,-1100]; 
P0 = point1_img_camera1;
P1 = [25,147,-800];
P2 = point1_img_camera2;
m = 1 ;
n=zeros(100,3);
n(m,:) = P0;
for m=1:100
    t = m/100;
    temp = (1-t)^2*P0+2*t*(1-t)*P1+t^2*P2;
    n = round(n);
    if (temp(1,1)-n(m,1))>=1
        n(m+1,:) = temp;
    end
    plot3(n(m,1),n(m,2),n(m,3),'o')
    hold on
end



clear
point1_img_camera1 = [25,147,-1000];  point2_img_camera1 = [238,147,-1100];  ...
        point1_img_camera2 = [100,147,-900];  point2_img_camera2 = [80,147,-1100]; 
P0 = point1_img_camera1;
P1 = [25,147,-800];
P2 = point1_img_camera2;
m = 1 ;
n(m,:) = P0;
while 1
    m = m+1;
    t = m/100;
    temp =(1-t)^2*P0+2*t*(1-t)*P1+t^2*P2;
    n = round(n);
    if temp(1,1) ~= n(m,1)
        n(m,:) = temp;
        plot3(n(m,1),n(m,2),n(m,3),'o')
        hold on
    end
    if temp(1,1) == P2(1,1)
        break
    end
end

clear
% 定义三个控制点
P0 = [25, -1000];
P1 = [25, -800];
P2 = [100, -900];

% 构建bezier曲线
t = linspace(0, 1, 1000);
x = (1-t).^2 * P0(1) + 2*(1-t).*t*P1(1) + t.^2*P2(1);
y = (1-t).^2 * P0(2) + 2*(1-t).*t*P1(2) + t.^2*P2(2);

% 提取x=25到x=100之间的点
idx = find(x>=25 & x<=100);
x_selected = round(x(idx)); % 四舍五入为整数
y_selected = round(y(idx));

% 绘制曲线
figure;
plot(P0(1), P0(2), 'o', 'MarkerSize', 10, 'LineWidth', 2); hold on;
plot(P1(1), P1(2), 'o', 'MarkerSize', 10, 'LineWidth', 2);
plot(P2(1), P2(2), 'o', 'MarkerSize', 10, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
plot(x_selected, y_selected, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
legend('P0', 'P1', 'P2', 'Bezier曲线', '选定点');

% 获取二值化图像的矩形边框
stats = regionprops('table', imfill_img_camera1, 'BoundingBox','Extrema');
a = stats.BoundingBox;
% 显示图像和边框
imshow(imfill_img_camera1);
hold on;
rectangle('Position', stats.BoundingBox, 'EdgeColor', 'r', 'LineWidth', 2);



figure(2)
plot3(point1_3d_img_camera1(1,1),point1_3d_img_camera1(1,2),point1_3d_img_camera1(1,3), 'bo', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot3(point2_3d_img_camera1(1,1),point2_3d_img_camera1(1,2),point2_3d_img_camera1(1,3), 'ro', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot3(point1_3d_img_camera2(1,1),point1_3d_img_camera2(1,2),point1_3d_img_camera2(1,3), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot3(point2_3d_img_camera2(1,1),point2_3d_img_camera2(1,2),point2_3d_img_camera2(1,3), 'go', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot3(P_a_camera1and2(1,1),P_a_camera1and2(1,2),P_a_camera1and2(1,3), 'y*', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot3(P_b_camera1and2(1,1),P_b_camera1and2(1,2),P_b_camera1and2(1,3), 'y*', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot3(P_c_camera1and2(1,1),P_c_camera1and2(1,2),P_c_camera1and2(1,3), 'y*', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot3(P_d_camera1and2(1,1),P_d_camera1and2(1,2),P_d_camera1and2(1,3), 'y*', 'MarkerSize', 10, 'LineWidth', 2)
axis equal
view(0,180)

    point1 = point1_3d_img_camera1(1,:)
    point2 = point1_3d_img_camera2(1,:)
    ab1 = tan(camera1)
    ab2 = -tan(camera2)
    cb1 = ab1*point1(1,1)+0*point1(1,3)
    cb2 = ab2*point2(1,1)+1*point2(1,3)
    line1 = [ab1, 0]
    line2 = [ab2, 1]
    Point = [line1;line2]\[cb1;cb2]
    Point_camera1and2 = [Point(1,1),point1(1,2),Point(2,1)]


A = magic(3);
B = [15; 15; 15];
x = A\B


