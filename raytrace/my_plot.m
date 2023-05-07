function my_plot(image,bubble,varargin)
% 本函数为气泡图像处理算法提供画圆、椭圆功能，使用方式如下：
% 
% 输入：
%     image：气泡图像
%     bubble：气泡边界信息
%     varargin：可变长度输入列表，用于指定画图内容
% 
% 输出：
%     无输出
% 
% 调用说明：
%     my_plot(image,bubble)
%                  默认画图方法，画所有圆和椭圆，不带框
%
% 
% 版本号VOL1.0，编写于2021年6月4日，作者：WG-Chen
%% default value


%% 计算部
figure;
imshow(image);
hold on;
[num2,~] = size(bubble);
for n=1:num2
    if(bubble{n,2}>1)
        [b,~]=size(bubble{n,6});
        for i=1:b
            if(bubble{n,6}{i,1}~=0)
                plot_ellipse(real(bubble{n,6}{i,1}(1,3)),real(bubble{n,6}{i,1}(1,4)),real(bubble{n,6}{i,1}(1,1)),real(bubble{n,6}{i,1}(1,2)),real(bubble{n,6}{i,1}(1,5)));
                hold on;
            end
        end
    elseif(bubble{n,2}==1 & bubble{n,6}~=0)
        plot_ellipse(real(bubble{n,6}(1,3)),real(bubble{n,6}(1,4)),real(bubble{n,6}(1,1)),real(bubble{n,6}(1,2)),real(bubble{n,6}(1,5)));
        hold on;
    elseif(bubble{n,2}==0.5)
        plot_ellipse(real(bubble{n,6}(2,2)),real(bubble{n,6}(2,1)),real(bubble{n,6}(1,1)),real(bubble{n,6}(1,1)),0);
        hold on;
    elseif(bubble{n,2}==0 & bubble{n,6}~=0)
        plot_ellipse(real(bubble{n,6}(1,3)),real(bubble{n,6}(1,4)),real(bubble{n,6}(1,1)),real(bubble{n,6}(1,2)),real(bubble{n,6}(1,5)));
        hold on;
    end
end
hold off

end


function plot_ellipse(x0,y0,a,b,theta)
%%  说明
%本程序画一个中心在（x0，y0)处的椭圆，其长短轴分别为a,b,椭圆沿Z轴转theta角
%%
num_t=1e2;
t=linspace(1,num_t+1,num_t);
x=x0+a*cos(2*pi/num_t*t)*cos(theta)-b*sin(2*pi/num_t*t)*sin(theta);
y=y0-a*cos(2*pi/num_t*t)*sin(theta)-b*sin(2*pi/num_t*t)*cos(theta);
global h1; %特殊用处，可删除
h1=plot(x,y,'r','LineWidth',1.5);  %特殊用处，可删除
end



