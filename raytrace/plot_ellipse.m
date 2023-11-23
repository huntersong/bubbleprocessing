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