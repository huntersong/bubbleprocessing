% 本函数的作用是为气泡图像处理算法提供基于霍夫变换的圆识别程序，使用方式如下：
% 
% 输入：
%     image：原始气泡图
%     varargin：可变长度输入列表，用于指定阈值阈值
% 
% 输出：
%     centers：识别到的圆心位置
%     raddis：识别到的半径
% 
% 调用说明：
%     [centers,raddis] = my_hough(image)
%                  默认预处理方式
%     [centers,raddis] = my_hough(image,varargin)
%                  指定阈值
% 
% 版本号VOL1.0，编写于2021年6月5日，作者：WG-Chen
function [centers,raddis] = my_hough(image,varargin)


% default value

image = adapthisteq(image);
[centers,raddis] = imfindcircles(image,[5 50],'ObjectPolarity','dark');


end